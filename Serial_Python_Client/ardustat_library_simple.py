import serial 
import socket
import glob
import pickle
from time import sleep,time
import subprocess
import sys
import os
import atexit

class ardustat:
	def __init__(self):	
		self.port = ""
		self.ser = serial.Serial()
		self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		self.mode = "serial"
		self.debug = False
		self.chatty = False
		self.groundvalue = 0
		self.p = None
		self.r_fixed_bool = False
		
	def kill_connection(self):
		try:
			print "Closing Time"
		except Exception as err:
			print "Latter!"
	
	def findPorts(self):
		"""A commands to find possible ardustat ports with no Arguments, """
		return glob.glob("/dev/tty.u*")
		
	def auto_connect(self,sport):
		"""This seems to work now!  Enter the port to try listening.  Try's to start a serial forwarder if one isn't started"""
		port = ""
		start_server = True
		#try:
		#	self.connect(sport)
		#except:
		#	start_server == True
			
		if start_server and os.name == "posix":
			#try os x
			if len(glob.glob("/dev/tty.u*")) > 0:
				port = glob.glob("/dev/tty.u*")
			elif len(glob.glob("/dev/ttyUSB*")) > 0:
				port = glob.glob("/dev/ttyUSB*")
			else:
				print "can't see any ardustats.  PEACE."
				sys.exit()
			if len(port) > 0:
				print port[0]
				try:
					self.p = subprocess.Popen(("python tcp_serial_redirect.py "+"-P "+str(sport)+" "+port[0]+" 57600").split())
					print "connected!"
				except:
					print "probably open trying to connect"
				sleep(3)
				self.s.connect(("localhost",sport))
				print "ardustat should be responding, trying a blink"
			else:
				print "not seeing any ardustats.  PEACE."
				sys.exit()
		
		
	def blink(self):
		self.rawwrite(" ");


	def connect(self,port=None):		
		if port == None:
			if len(glob.glob("/dev/tty.u*")) > 0:
				port = glob.glob("/dev/tty.u*")
			elif len(glob.glob("/dev/ttyUSB*")) > 0:
				port = glob.glob("/dev/ttyUSB*")
			else:
				print "can't see any ardustats.  PEACE."
				sys.exit()
		port = port[0]
		if self.mode == "serial":
			self.ser = serial.Serial(port,57600)
			self.ser.timeout = .1
			self.ser.open()
			return "connected to serial"
		if self.mode == "socket":
			self.s.connect(("localhost",port))
			return "connected to socket"
	
	def ocv(self):
		self.rawwrite("-0000")
	
	def potentiostat(self,potential):
		"""Argument: Potential (V).  Sets the potentiostat"""
		#potential = potential#+self.groundvalue
		if potential < 0: potential = str(2000+int(1023*(abs(potential)/5.0))).rjust(4,"0")
		else: potential = str(int(1023*(potential/5.0))).rjust(4,"0")
		if potential == "2000": potential = "0000"
		self.rawwrite("p"+potential)
	
		
	def rawwrite(self,command):
		if self.mode == "serial":
			self.ser.write(command)
		if self.mode == "socket":
			self.s.send(command)
			sleep(.05)
			if self.chatty:
				return self.parsedread()
				
	
	def rawread(self):
		self.blink()
		self.ser.flushInput()
		sleep(.01)
		self.rawwrite("s0000")
		self.rawwrite("s0000")
		self.rawwrite("s0000")
		sleep(.01)
		if self.mode == "serial":
			a = self.ser.read(1)			
			while 1:
				a += self.ser.read(1)
				#print a
				if a.find("ST\r\n") > 1: return a.strip()
			
		if self.mode == "socket":
			a = ""
			while 1:
				a += self.s.recv(1024)
				if a.find("ST\r\n") > 1: return a.strip()

	def solve_for_r(self,input_r,v_in,v_out):
		return input_r*((v_in/v_out)-1)
	
	
	def galvanostat(self,current):
		"""Tries to pick the ideal resistance and sets a current difference"""
		#V = I R -> I = delta V / R
		#goal -> delta V = .2 V
		R_goal = abs(.1 / current)
		R_real = 10000
		R_set = 0
		err = 1000
		for d in self.res_table:
			this_err = abs(self.res_table[d][0]-R_goal)
			if this_err < err:
				err = this_err
				R_set = d
				R_real = self.res_table[d][0]
		#Solve for real delta V
		if self.r_fixed_bool:
			R_real = r_fixed
		delta_V = abs(current*R_real)
		if self.debug: print current,delta_V,R_real,R_set
		potential = str(int(1023*(delta_V/5.0))).rjust(4,"0")
		if current < 0:
			potential = str(int(potential)+2000)
		if self.debug: print "gstat setting:", potential
		print potential
		if potential == "2000" or potential == "0000":
			self.ocv()
		else:
			#Write!
			self.rawwrite("r"+str(R_set).rjust(4,"0"))
			sleep(.1)
			self.rawwrite("g"+str(potential))
		
		

	def load_resistance_table(self,id):
		self.res_table =  pickle.load(open("unit_"+str(id)+".pickle"))
		if self.debug: print self.res_table
		self.id = id

	def calibrate(self,known_r,id):
		"""Runs a calibration by setting the resistance against a known resistor and solving the voltage divider equation.  Cycles through all possible resistances"""
		
		#Make a Big List of Correlations
		ressers = []
		self.rawwrite("R")
		sleep(.1)
		self.rawwrite("r0000")		
		for i in range(0,10):
			for y in range(0,255):
				self.rawwrite("r"+str(y).rjust(4,"0"))
				sleep(.05)
				values = self.parsedread()
				if values['valid']:
					solved_r = self.solve_for_r(known_r,values['DAC0_ADC'],values['cell_ADC'])
					print values['pot_step'], solved_r
					ressers.append([int(values['pot_step']), solved_r])
				else: print "bad read"
					
		self.ocv()
		
		#Make a Big List of Correlations
		big_dict = {}
		for r in ressers:
			try:
				big_dict[r[0]].append(r[1])
			except:
				big_dict[r[0]] = []
				big_dict[r[0]].append(r[1])
				
		#Find Values
		final_dict = {}
		for b in big_dict.keys():
			final_dict[b] = [sum(big_dict[b])/len(big_dict[b]),(max(big_dict[b])-min(big_dict[b]))/2.0]
		pickle.dump(final_dict,open("unit_"+str(id)+".pickle","w"))
		return final_dict
	
	def parsedread(self):
		return self.parseline(self.rawread())
		
	def last_raw_reading(self):
		out = ""
		while self.ser.inWaiting() > 0: out = self.ser.readline()
		return out
	
	
	def moveground(self):
		"""Argument: Potential (V).  Moves ground to allow negative potentials.  Check jumper position!"""
		potential = str(int(1023*(self.groundvalue/5.0))).rjust(4,"0")
		self.rawwrite("d"+potential)
	
	def refbasis(self,reading,ref):
		"""Argument: raw ADC reading, raw ADC basis.  Returns an absolute potential based on the ADC reading against the 2.5 V reference
		   (reading from pot as a value between 0 and 1023, reference value in V (e.g. 2.5))"""
		return round((float(reading)/float(ref))*2.5,3)
	
	def resbasis(self,reading,pot):
		"""Argument, raw pot setting, max pOT reading.   Returns the value for the givening potentiometer setting 
			(reading as value between 0 and 255, pot lookup variable).  Wildly Inaccurate.  Don't use."""
		return round(10+(float(reading)/255.0)*pot,2)
	
	def parseline(self,reading):
		outdict = {}
		#format GO,1023,88,88,255,0,1,-0000,0,0,510,ST 
			#splitline[0] -> "GO"
			#splitline[1] -> DAC Setting
			#splitline[2] -> Cell Voltage
			#splitline[3] -> DAC Measurement
			#splitline[4] -> DVR Setting
			#splitline[5] -> "setout" variable in firmware. Don't know what this represents
			#splitline[6] -> Mode (0 = manual, 1 = OCV, 2 = potentiostat, 3 = galvanostat)
			#splitline[7] -> Last command received
			#splitline[8] -> GND Measurement
			#splitline[9] -> Reference Electrode
			#splitline[10] -> Reference Voltage
			#splitline[11] -> "ST"
		try:
			if reading.find("GO") == 0 and reading.find("ST") and reading.rfind("GO") == 0:
				outdict['valid'] = True
				outdict['raw'] = reading
				outdict['time'] = time()
				parts = reading.split(",")
				outdict['ref'] = float(parts[len(parts)-3])
				outdict['DAC0_ADC'] = self.refbasis(parts[3],outdict['ref'])-self.refbasis(parts[8],outdict['ref'])
				outdict['cell_ADC'] = self.refbasis(parts[2],outdict['ref'])-self.refbasis(parts[8],outdict['ref'])
				outdict['ref_ADC'] = self.refbasis(parts[9],outdict['ref'])-self.refbasis(parts[8],outdict['ref'])
				outdict['pot_step'] = parts[4]
				outdict['work_v_ref'] = outdict['cell_ADC'] - outdict['ref_ADC']
				outdict['unit_no'] = int(parts[11])
				
				##Try to read from the res_table, otherwise make the dangerous assumption
				try:
					outdict['res'] = self.res_table[int(outdict['pot_step'])][0]
				except Exception, err:
					if self.debug: print err
					outdict['res'] = self.resbasis(outdict['pot_step'],10000.0)
				outdict['current'] = (float(outdict['DAC0_ADC'])-float(outdict['cell_ADC']))/outdict['res']			
			else:
				outdict['valid'] = False
				outdict['raw'] = reading
				
		except Exception as error:
			if self.debug: print error
			outdict['valid'] = False
			outdict['raw'] = reading
			
		return outdict
