import numpy
import random

class ExponentialPlugin:
   def input(self, filename):
      filestuff = open(filename, 'r')
      firstline = filestuff.readline()
      self.bacteria = firstline.split(',')
      if (self.bacteria.count('\"\"') != 0):
         self.bacteria.remove('\"\"')
      self.n = len(self.bacteria)

      self.ADJ = numpy.matrix(numpy.zeros([self.n, self.n]))
      self.P = numpy.matrix(numpy.zeros([self.n, 1]))
 
      i = 0
      for line in filestuff:
         contents = line.split(',')
         for j in range(self.n):
            if (i != j): 
              self.ADJ[i,j] = float(contents[j+1])
         i += 1

      for i in range(self.n):
         self.bacteria[i] = self.bacteria[i].strip()
         self.P[i,0] = 1.0/self.n
     

      self.EPS = random.random()*0.000001
      self.MU = (numpy.amax(self.ADJ) - numpy.amin(self.ADJ)) / 2.0 + self.EPS 


   def run(self):
      oldP = numpy.zeros([self.n, 1])
      iteration = 0
      while (numpy.amax(self.P - oldP) > self.EPS):
         oldP = self.P.copy()
         tmp = numpy.exp(1.0/self.MU)*numpy.transpose(self.ADJ)*self.P
         self.P = tmp / numpy.linalg.norm(tmp)
         iteration += 1
      self.k = numpy.transpose(self.ADJ)*self.P
      print "*** Converged after ", iteration, " iterations ***"

   def output(self, filename):
     UG = []
     for i in range(self.n):
      UG.append((abs(self.k[i,0]), self.bacteria[i]))
     UG.sort()
     UG.reverse()
     # Formatted for Cytoscape, borrowed from ATria
     outfile = open(filename, 'w')
     outfile.write("Name\tCentrality\tRank\n")
     centvals = numpy.zeros([len(UG)])
     for i in range(len(UG)):
       print (UG[i][1], UG[i][0])
       bac = UG[i][1]
       if (bac[0] == '\"'):
          bac = bac[1:len(bac)-1]
       if (UG[i][0] != UG[len(UG)-1][0]):
         outfile.write(bac+"\t"+str(abs(UG[i][0]))+"\t"+str(len(UG)-i)+"\n")
       else:
         outfile.write(bac+"\t"+str(abs(UG[i][0]))+"\t"+"0\n")
       centvals[i] = abs(UG[i][0])

     print "Wrote file: ", file
     print "Min centrality: ", numpy.min(centvals)
     print "Max centrality: ", numpy.max(centvals)
     mymean = numpy.mean(centvals)
     stddev = numpy.std(centvals)
     print "Standard Deviation: ", stddev
     print "Two STDs back: ", mymean - 2*stddev
     print "One STD back: ", mymean - stddev
     print "One STD forward: ", mymean + stddev
     print "Two STDs forward: ", mymean + 2*stddev

