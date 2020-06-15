
print '?!'*21
print 'Welcome to the sBraid production line'
print 'Nathan Borggren'
print 'August 2009'
print '?!'*21

z = raw_input("How many magnetic lines do you want?")
mlines = int(z)



print " "
print "The Rules of the Game:"
print "alphabet:"
print "    \'mn\' crosses line n, n+1"
print "    \'Mn\' is the inverse of mn, crossing line n+1 over n"
print "    \'xn'\' makes a copy of line n and switches line n+1 to n XOR n+1" 
print "    \'Xn'\' makes a copy of line n+1 and switches line n to n XOR n+1" 
print " "


eBraid = raw_input("how should I braid them?")
