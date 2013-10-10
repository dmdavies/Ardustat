from Ardustat_Class import Everything as Node
#-------------------------------------------------------------------------#           

#setup 
#-------------------------------------------------------------------------#
hostsite = 'http://localhost:4000/'
n = Node(hostsite)
ardustat_id = 24
file_name = "Galvo_test_:)_"

#set parameters
#-------------------------------------------------------------------------#
cycles = 1
DAC2_value = 1
min_current = -.001 #amps
max_current = .001 #amps
rate = 2 #mV/s
read_delay = .5 #second
time_delay = 2 #second
increment_size = .00007 #amps
debug = False
print "greg is a bogan aah"

n.galvanostat_run(file_name,ardustat_id,min_current,max_current,increment_size,time_delay,read_delay,cycles,DAC2_value)
       

