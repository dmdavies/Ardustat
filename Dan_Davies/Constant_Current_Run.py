from Ardustat_Class import Everything as Node
#-------------------------------------------------------------------------#           

#setup 
#-------------------------------------------------------------------------#
hostsite = 'http://localhost:3000/'
n = Node(hostsite)
ardustat_id = 22
file_name = "Constant_current_test_"

#set parameters
#-------------------------------------------------------------------------#
cycles = 1
DAC2_value = 1
min_voltage = -10 #V
max_voltage = 10 #V
constant_current = .0005 #Amps make sure this is a positive value - python will do discharge and charge shit
read_delay = .5 #second
capacity_limit =100 #Farads
print "start"
       
#new run :)
#write now, this only calculates a very approximate capacity. will write a fancier one
#that calculates it more accurately..

n.constant_current_cycle_voltage(file_name,ardustat_id,min_voltage,max_voltage,constant_current,read_delay,cycles,DAC2_value,capacity_limit)
print file_name
print "done."


