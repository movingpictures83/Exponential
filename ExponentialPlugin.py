import numpy
import random
random.seed(1234)   # ONLY FOR TESTING, COMMENT OUT FOR BETTER RANDOMNESS
import PyPluMA

EPS=1e-8

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
      PyPluMA.log("*** Converged after "+str(iteration)+" iterations ***")

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
       PyPluMA.log(str((UG[i][1], UG[i][0])))
       bac = UG[i][1]
       rank = len(UG)-i
       if (i != 0 and abs(UG[i][0]-UG[i-1][0]) < EPS):
          rank = oldrank
       oldrank = rank
       if (bac[0] == '\"'):
          bac = bac[1:len(bac)-1]
       #if (UG[i][0] != UG[len(UG)-1][0]):
       outfile.write(bac+"\t"+str(abs(UG[i][0]))+"\t"+str(rank)+"\n")
       #else:
       #  outfile.write(bac+"\t"+str(abs(UG[i][0]))+"\t"+"0\n")
       centvals[i] = abs(UG[i][0])

     PyPluMA.log("Wrote file: "+filename)
     PyPluMA.log("Min centrality: "+str(numpy.min(centvals)))
     PyPluMA.log("Max centrality: "+str(numpy.max(centvals)))
     mymean = numpy.mean(centvals)
     stddev = numpy.std(centvals)
     PyPluMA.log("Standard Deviation: "+str(stddev))
     PyPluMA.log("Two STDs back: "+str(mymean - 2*stddev))
     PyPluMA.log("One STD back: "+str(mymean - stddev))
     PyPluMA.log("One STD forward: "+str(mymean + stddev))
     PyPluMA.log("Two STDs forward: "+str(mymean + 2*stddev))

