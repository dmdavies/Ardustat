import urllib2
import pickle
import time
from time import sleep
import pylab as pl
import numpy as np

#made things negative - need to do relative things for stuff other than CV.


class Everything():
    def __init__(self,base):
        self.base = base
        self.bwrite = self.base+"write/"
        self.bread = self.base+"read/"
        self.debug = False
        
    def sread(self):
        foo = urllib2.urlopen(self.bwrite+"s")
        sleep(.05)
        foo = urllib2.urlopen(self.bread)
        buff = foo.read()
        out = buff[buff.rfind("GO"):buff.rfind("ST")+2]
        return out
    
    def swrite(self,stringer):
        foo = urllib2.urlopen(self.bwrite+stringer)
        sleep(.05)
        return foo.read()
    
    #clean this up
    def parseline(self,reading,res_table):
        outdict = {}
        #print "reading: ",reading
        outdict['valid'] = True
        parts = reading.split(",")
        outdict['ref'] =float(parts[-2])
        outdict['DAC0_ADC'] = self.refbasis(parts[3],outdict['ref'])
        outdict['cell_ADC'] = self.refbasis(parts[2],outdict['ref'])
        outdict['pot_step'] = parts[4]
        outdict['other_DAC'] = self.refbasis(parts[-4],outdict['ref'])
        outdict['ref_electrode'] = self.refbasis(parts[-3],outdict['ref'])
        #print 'ref_electrode: ',outdict['ref_electrode']
       
        #making potentials relative to DAC2
        outdict['tot_potential'] = outdict['DAC0_ADC'] - outdict['other_DAC'] #making potentials relative to DAC2
        outdict['cell_potential'] = outdict['cell_ADC'] - outdict['other_DAC'] #making potentials relative to DAC2
        outdict['res'] = res_table[int(outdict['pot_step'])][0]
        outdict['current'] = (float(outdict['tot_potential'])-float(outdict['cell_potential']))/outdict['res']
        outdict['pot-ref'] = (float(outdict['cell_potential'])-float(outdict['ref_electrode']))
        outdict['adc-ref'] = (float(outdict['cell_ADC'])-float(outdict['ref_electrode']))
        if self.debug:
            print outdict
        return outdict
        
    def load_resistance_table(self,ardustat_id): 
        #don't need the self. here i dont think...?? need to put in an id
        if (ardustat_id == 22):
            self.res_table = pickle.load(open("res_table_full_22.p")) #res_table_full refers to the res_table for ardustat 24
        elif (ardustat_id == 24):
            self.res_table = pickle.load(open("res_table_full.p"))
            self.res_table[255] = [10165,5] 
        else: print "error - wrong ardustat id"
        return self.res_table
        
    def refbasis(self,reading,ref):
        return round((float(reading)/float(ref))*2.5,3)
        
    def ocv(self):
        self.swrite("-0000")

    def potentiostat(self,potential):
        voltage = str(int(1023*(abs(potential)/5.0))).rjust(4,"0")
        if potential < 0:
            voltage = int(voltage) + 2000
            voltage = str(voltage)
        self.swrite("p"+voltage) 
        
    def parsedread(self,res_table):
        return self.parseline(self.sread(),res_table)
    
       #when galvonstat is called - to achieve a negative current a negative
       #potential must be called. to do this - is this coded into the firmware?
       #i.e. when gXXXX is called - does it understand to set the potential to 
       #below the other_DAC...   
    def galvanostat(self,current,res_table):
        #print "trying some things"
        """Tries to pick the ideal resistance and sets a current difference"""
        #V = I R -> I = delta V / R
  		#goal -> delta V = .2 V
        R_goal = .1 / current
        R_real = 10000
        R_set = 0
        err = 1000
        for d in res_table:
       	    this_err = abs(res_table[d][0]-R_goal)
       	    if this_err < err:
      		    err = this_err
      		    R_set = d
      		    R_real = res_table[d][0]
        #Solve for real delta V
        delta_V = abs(current*R_real)
        #print "delta_V:",delta_V, "R_real:",R_real
        potential = str(int(1023*(delta_V/5.0))).rjust(4,"0")
        if current < 0:
       	    potential = str(int(potential)+2000)
        print "gstat setting:", potential
        print "rstat setting:", R_set
        self.swrite("r"+str(R_set).rjust(4,"0"))
        #s.readlines()
        time.sleep(.1)
        self.swrite("g"+str(potential)) 
        #s.readlines()
        print "resistance:",R_real

    
    def appender_test(self,reading): #need to modify this so that it changes the file you write to automatically
        print reading['cell_potential'],reading['current']
        open('C:/Python27/projects/ardustat/document.csv','a').write(str(reading['cell_ADC'])+',' + str(reading['current']) + '\n')
    
    #need to change to make better
    def galv_appender(self,reading,cycle,file_name,time_start):
        rawreading = self.sread()
        print rawreading
        print 'cell potential: ', reading['cell_potential']
       	print 'cell_pot - ref_pot: %r' %reading['pot-ref'], 'current: %r' %reading['current']
       	tdiff = str(time.time()-time_start)
       	out = tdiff+","+str(reading['pot-ref'])+","+str(reading['current'])+","+str(reading['cell_potential'])+str(cycle)+"\n"
       	open(file_name,"a").write(out)
       	
    def constant_current_appender(self,file_name,cycle,time_start,voltage,capacity,current):
        print 'start appender'
        rawreading = self.sread()
        print rawreading
        print "current: ",current
        print "voltage: ",voltage
        print "capacity: ",capacity
        print "appender donezo"
        tdiff = str(time.time()-time_start)
       	out = tdiff+","+str(current)+","+str(voltage)+","+str(capacity)+str(cycle)+"\n"
       	open(file_name,"a").write(out)
       	#write a plotting function for this too.
        
       	
    def constant_current_cycle_voltage(self,file_name,ardustat_id,min_voltage,max_voltage,constant_current,read_delay,cycles,DAC2_value,capacity_limit):
        time_start = time.time()
        file_name = file_name+"_time"+str(int(time_start))+".dat"
        res_table = self.load_resistance_table(ardustat_id)
        self.set_DAC2(DAC2_value)
        read = self.parsedread(res_table)
        print "dac2: %r" %read['other_DAC']
        cycle = 0
        
        #discharge
        while cycle < cycles:
            capacity = 0
            print 'cycle: ' ,cycle
            voltage = max_voltage
            discharge_current = -1 * constant_current
            time_step = time.time()
            while (voltage > min_voltage) or (capacity < capacity_limit):
                self.galvanostat(discharge_current,res_table)
                print "attempted current: ",discharge_current
                time.sleep(read_delay)
                read = self.parsedread(res_table)
                voltage = read['cell_potential'] #i think
                time_dif = time.time() - time_step
                print "time difference: ",time_dif
                capacity = read['current'] * time_dif / voltage#current * time / volatge 
                current = read['current']
                self.constant_current_appender(file_name,cycle,time_start,voltage,capacity,current)
                
        #charge 
            capacity = 0
            time_step = time.time()
            voltage = min_voltage
            while voltage < max_voltage or capacity < capacity_limit:
                print "attempted current: ", constant_current
                self.galvanostat(constant_current,res_table)
                time.sleep(read_delay)
                read = self.parsedread(res_table)
                voltage = read['cell_potential'] #i think
                time_dif = time.time() - time_step
                capacity = read['current'] * time_dif / voltage#current * time / volatge 
                current = read['current']
                self.constant_current_appender(file_name,cycle,time_start,voltage,capacity,current)
                
            cycle = cycle +1
            
        return file_name
                

            
        
        
    def galvanostat_run_time_based(self,file_name,ardustat_id,min_current,max_current,increment_size,time_delay,read_delay,cycles,DAC2_value):
        #should maybe write a bunch of if loops for people to select what they want to measure.
        #write the plotting functions based on the settings people gave us..
        
        amps = min_current
        cycle = 0
        res_table = self.load_resistance_table(ardustat_id)
        self.set_DAC2(DAC2_value)
        read = self.parsedread(res_table)
        print "dac2: %r" %read['other_DAC']
        time_start = time.time()
        
        
        #not too sure what is going on with the order and mixing up of attempt potential and the potential reading - will need to sort out
        while cycle < cycles:
       	#Scan Up
       	    
       	    print 'cycle:',cycle
            amps = min_current
       	    while amps < max_current:
       	        print "attempt current: %r" %amps
       	        step_time = time.time()
      		self.galvanostat(amps,res_table)
      		amps = amps + increment_size
      		a = 0
      		while (a-step_time) < time_delay:
     			time.sleep(read_delay)
     			read = self.parsedread(res_table)
     			self.galv_appender(read,cycle,file_name,time_start)
     			a = time.time()
     	#scan down
     	    amps = max_current
     	    while amps > min_current:
       	        print "attempt current: %r" %amps
       	        step_time = time.time()
      		self.galvanostat(amps,res_table)
      		amps = amps - increment_size
      		a = 0
      		while (a-step_time) < time_delay:
     			time.sleep(read_delay)
     			read = self.parsedread(res_table)
     			self.galv_appender(read,cycle,file_name,time_start)
     			a = time.time()
            cycle = cycle +1 
       	return file_name #return t
       	
    def calibrate_test(self,known_res):
        self.swrite("R")
        self.swrite("r0255")
        values = self.parsedread()
        print 'values'
        print values
        built_in_res = self.get_r(known_res,values['DAC0_ADC'],values['cell_ADC'])
        print '%r ohms' % built_in_res

    def get_r(self,known_res,DAC,ADC):
        return float(known_res)*((float(DAC)/float(ADC))-1)
        
    
    def calibrate(self,known_res):
        
        self.swrite("R")
        self.swrite("r0001")
        ressers = []
        values = self.parsedread()
        print ' here we go '
        for i in range(0,10):
            print i 
            for y in range(0,256):
                self.swrite("r"+str(y).rjust(4,"0"))
                values = self.parsedread()
                built_in_res = self.get_r(known_res,values['DAC0_ADC'],values['cell_ADC'])
                #print values['pot_step'],built_in_res
                open('C:/Python27/projects/ardustat/bogan8.csv','a').write(str(values['pot_step'])+',' + str(built_in_res) + '\n') 
                ressers.append([int(values['pot_step']), built_in_res])    
        #values = parsedread()
        
        #Make a Big List of Correlations
  		
        big_dict = {}
        for r in ressers:
   	    try:
    	       big_dict[r[0]].append(r[1])
   	    except:
    	       big_dict[r[0]] = []
    	       big_dict[r[0]].append(r[1])
    	   
        print big_dict
    				
        #Find Values
        final_dict = {}
        print len(big_dict.keys())
        print big_dict.keys()
        for b in big_dict.keys():
            print b
            final_dict[b] = [sum(big_dict[b])/len(big_dict[b]),(max(big_dict[b])-min(big_dict[b]))/2.0]
   	pickle.dump(final_dict,open("res_table_full_22.p","wb"))
        return final_dict
        
    def CV_potentiostat(self,values,potential):
        #can surely write this so that it is faster and neater but i think it works.
        print "potential trying to reach: ",potential
        voltage = str(int(1023*(abs(potential)/5.0))).rjust(4,"0")
        if potential < 0:
            voltage = int(voltage) + 2000
            voltage = str(voltage)
        #print "what im actually sending it: p",DAC1
        print 'sending arduino: ',("p"+voltage)
        self.swrite("p"+voltage)
        
    def cyclic_voltammagram_time_test(self,cycles,min_potential,max_potential,rate,ardustat_id,file_name,read_delay,DAC2_value):
        
        res_table = self.load_resistance_table(ardustat_id)
        
        time_start = time.time()
        cycle = 0
        #self.CV_sleep()
        file_name_raw = file_name+"_time"+str(int(time_start))+"_raw.dat"
        file_name = file_name+"_time"+str(int(time_start))+".dat"
        
        
        self.set_DAC2(DAC2_value)
        reading = self.parsedread(res_table)
        print "dac2: %r" %reading['other_DAC']
        
        #not too sure what is going on with the order and mixing up of attempt potential and the potential reading - will need to sort out
        while cycle < cycles:
       	#Scan Up
       	    print 'cycle:',cycle
            output = min_potential
            i = 1
       	    while i < 20:
       	        print "attempt potential: %r" %output
       	        step_time = time.time()
      		self.swrite('p512')
      		i = i+1
      		while (time.time()-step_time) < (5/rate):
     			time.sleep(read_delay)
     			read = self.parsedread(res_table)
     			self.CV_appender_time(read,cycle,file_name,file_name_raw,time_start)
     			
        
    def cyclic_voltammagram_time(self,cycles,min_potential,max_potential,rate,ardustat_id,file_name,read_delay,DAC2_value):
        
        res_table = self.load_resistance_table(ardustat_id)
        
        time_start = time.time()
        cycle = 0
        #self.CV_sleep()
        file_name_raw = file_name+"_time"+str(int(time_start))+"_raw.dat"
        file_name = file_name+"_time"+str(int(time_start))+".dat"
        
        #will change to allow user to set DAC2 themselves - for now they are stupid though.
        self.set_DAC2(DAC2_value)
        read = self.parsedread(res_table)
        #print "dac2: %r" %read['other_DAC']
        
        #not too sure what is going on with the order and mixing up of attempt potential and the potential reading - will need to sort out
        while cycle < cycles:
       	#Scan Up
       	    #print 'cycle:',cycle
            output = min_potential
       	    while output < max_potential:
       	        print "attempt potential: %r" %output
       	        step_time = time.time()
      		self.CV_potentiostat(read,output)
      		output = output
      		while (time.time()-step_time) < (5/rate):
     			time.sleep(read_delay)
     			read = self.parsedread(res_table)
     			self.CV_appender_time(read,cycle,file_name,file_name_raw,time_start)
     			
     	    output = max_potential
     	    while output > min_potential:
       	        print "attempt potential: %r" %output
       	        step_time = time.time()
      		self.CV_potentiostat(read,output)
      		output = output - .005
      		while (time.time()-step_time) < (5/rate):
     			time.sleep(read_delay)
     			read = self.parsedread(res_table)
     			self.CV_appender_time(read,cycle,file_name,file_name_raw,time_start)
            cycle = cycle +1 
       	return file_name   #return the file_name for easy-peasy printing :)
       	    
       	
#this one seems dodgy - not too worried cause pretty sure we aren't going to use it anyway.. 
#fix this to make it like run properly like the time one...
     
    def cyclic_voltammagram_step(self,cycles,min_potential,max_potential,ardustat_id,file_name):
        
        res_table = self.load_resistance_table(ardustat_id)
        #create arrays + a function for logging data
        cycle = 0
        file_name = file_name+"_"+"step"+".dat"
        #self.CV_sleep(res_table,file_name)
        
        while cycle < cycles:
            i = min_potential
       	    while i < max_potential: 
       	        print "potential: %r" %i 
       	        self.potentiostat(i)
       	        time.sleep(.5)
       	        reading = self.parsedread(res_table)
       	        self.CV_appender_step(reading,file_name) 
       	        i = i + .05  
       	    
       	    i = max_potential
       	    while i > min_potential: 
       	        print "potential: %r" %i 
       	        self.potentiostat(i)
       	        time.sleep(.5)
       	        reading = self.parsedread(res_table)
       	        self.CV_appender_step(reading,file_name) 
       	        i = i - .05    
       	    cycle = cycle +1 
        return file_name   #return the file_name for easy-peasy printing :)
        
    def CV_appender_time(self,reading,cycle,file_name,file_name_raw,time_start):
        print 'start'
        rawreading = self.sread()
        print rawreading
        print 'cell potential: ', reading['cell_potential']
        print 'cell potential - reference potential: %r' %reading['pot-ref']
       	print 'current: %r' %reading['current']
       	print 'ref-electrode',reading['ref_electrode']
       	print 'adc measurement: ',reading['cell_ADC']
       	print 'adc-ref',reading['adc-ref'] 
       	print '\n'
       	tdiff = str(time.time()-time_start)
       	out = tdiff+","+str(reading['adc-ref'])+","+str(reading['current'])+","+str(reading['cell_potential'])+str(cycle)+"\n"
       	#out2= rawreading+"\n"
       	open(file_name,"a").write(out)
       	#open(file_name_raw,"a").write(out2)
       	
       	
    def CV_appender_step(self,reading,file_name):
        print 'cell potential: %r' %reading['cell_potential'], 'current: %r' %reading['current']
        open(file_name,'a').write(str(reading['cell_potenial'])+',' + str(reading['current']) + '\n')
        
    def CV_sleep(self): #Allows cell to settle and picks starting potential based on OCV
        self.ocv()
        for i in range(0,60):
       	    time.sleep(1)
       	
       	
    def plot_CV_time(self,file_name):
        print "kevin smells like he and brooker had a good time last night"
        
        foo = np.genfromtxt(file_name,skip_header=0,delimiter=",").transpose()
        print foo[1]
        print foo[2]*1000 #convert from amps to milliamps
        pl.plot(foo[1],foo[2]*1000)    
        pl.legend()
        pl.title(file_name)
        pl.xlabel('voltage, V')
        pl.ylabel('current, mA')
        pl.show()
        
    def plot_CV_fixer(self,file_name):
        foo = np.genfromtxt(file_name,skip_header=0,delimiter=",").transpose()
        print foo[1]
        print foo[2]*1000 #convert from amps to milliamps
        
        
        

        #find a way to save the plot to file.
    def plot_CV_step(self,file_name):
        print "start"
        
        #datalist = pl.loadtxt(file_name)
        #for data in datalist:
        #    pl.plot(data[:,1], data[:,2], '-')
        foo = np.genfromtxt(file_name,skip_header=0,delimiter=",").transpose()
        print foo[0]
        print foo[1]
        pl.plot(foo[0],foo[1])    
        pl.legend()
        pl.title(file_name+'pot-ref')
        pl.xlabel('voltage')
        pl.ylabel('current')
        pl.show()
        
        
        
    def set_DAC2(self,potential):
        potential = str(int(1023*(potential/5.0))).rjust(4,"0")
        self.swrite("d"+potential)
        
        
        
    
         
        
