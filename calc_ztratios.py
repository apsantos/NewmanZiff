#!/usr/bin/env python
"""moduli

 Calculate the Ziff-Torquato ratios, R1 and R2:
    Ziff, R. M., & Torquato, S. (2017). Percolation of disordered jammed sphere packings. J. Phys. A: Math. Theor., 50, 0850011. https://doi.org/10.1088/1751-8121/aa5664
 							                                    
 -Andrew P. Santos					                            
 
"""
import sys, argparse
import numpy as np

class NewmanZiff(object):
    """Callable moduli
    """
    def __init__(self):
        self.p = []
        self.big = []
        self.big2 = []
        self.moment = []

    def addParser(self, parser):
        """
        Get relevant values from an argparse parser
        """
        if (parser.parse_args().start_line != None):
            self.start_line = parser.parse_args().start_line
        else:
            print 'assuming startline is 1'
            self.start_line = 1

        if parser.parse_args().filename:
            self.filenames = []
            self.filenames.append( parser.parse_args().filename )
        elif parser.parse_args().filenamelist:
            self.filenames = parser.parse_args().filenamelist
            #self.filenames = self.readFilelist(parser.parse_args().filenamelist)
        else:
            return -1
            
        if parser.parse_args().outfilename:
            self.ofilename = parser.parse_args().outfilename

        if parser.parse_args().occupationfraction_collumn >= 0:
            self.ic_p = parser.parse_args().occupationfraction_collumn

        if parser.parse_args().biggestcluster_collumn >= 0:
            self.ic_big = parser.parse_args().biggestcluster_collumn

        if parser.parse_args().biggestcluster2_collumn >= 0:
            self.ic_big2 = parser.parse_args().biggestcluster2_collumn

        if parser.parse_args().moment_collumn >= 0:
            self.ic_moment = parser.parse_args().moment_collumn

    def read(self, filelist=None):
        if filelist == None:
            filelist = self.filenames

        for ifilename in filelist:
            self.readOccupation( ifilename )
            self.convolution()
            self.calcRatios()

    def readOccupation(self, filename=None):
        """
        Read the strain and stress
        """
        if filename == None:
            filename = self.filenames[0]

        try:
            ifile = open(filename, 'r')
        except IOError:
            raise IOError('cannot find: %s. Either put the file in the working directory or fix the input file.' % filename)

        max_ic = max( self.ic_p, self.ic_big, self.ic_big2, self.ic_moment)
        i_line = 0
        for line in ifile:
            i_line += 1
            if (i_line < self.start_line):
                continue

            data = line.strip().split()
            if len(data) <= max_ic:
                raise IOError('Not enough collumns in %s based on input information.' % filename)

            self.p.append( float( data[self.ic_p] ) )
            self.big.append( float( data[self.ic_big] ) )
            self.big2.append( float( data[self.ic_big2] ) )
            self.moment.append( float( data[self.ic_moment] ) )

    def convolution(self):
        N = len(self.p)
        B = np.ones( (N, N, N) )
        for i in range( N ):
            ip = self.p[i]
            nmax = ip*N
            for n in range( N ):
                if ( n > nmax):
                    B[N,n,i] = B[N,n-1,i] * (N - n + 1) * ip / float(n * (1 - ip))
                elif ( n < nmax):
                    B[N,n,i] = B[N,n+1,i] * (n + 1) * (1 - ip) / float((N - n) * ip)


    def calcRatios(self):
        from scipy import stats
        # get the slope
        max_strain = max(self.strain)
        min_strain = min(self.strain)
        tmp_strain = []
        tmp_stress = []
        ndata = len(self.strain)
        for i in range(ndata):
            if self.find_window:
                tmp_strain.append( self.strain[i] )
                tmp_stress.append( self.stress[i] )
            else:
                if ( (self.strain[i] < self.window) or (self.strain[i] > -self.window)):
                    tmp_strain.append( self.strain[i] )
                    tmp_stress.append( self.stress[i] )

        if self.find_window:
            niteration = 100
            # start with the full range of strain data, and parse down until the data is 'linear' enough
            for i in range(niteration):
                slope, intercept, R, p, std = stats.linregress( np.array(tmp_strain), np.array(tmp_stress) ) 
                if R**2.0 > 0.9995:
                    self.moduli.append(slope)
                    self.std_moduli.append(std)
                    return
                else:
                    max_strain *= 0.9
                    min_strain *= 0.9
                    tmp_strain = []
                    tmp_stress = []
                    for i in range(ndata):
                        if ( (self.strain[i] < max_strain) and (self.strain[i] > min_strain)):
                            tmp_strain.append( self.strain[i] )
                            tmp_stress.append( self.stress[i] )
        else:
            [intercept, slope], R2 = np.polyfit( np.array(tmp_strain), np.array(tmp_stress), 1)

    def write(self, filename=None):
        """
        Write the arrest times
        """
        if filename == None:
            filename = self.ofilename

        ofile = open(filename, 'w')

        ofile.write('# file moduli confidence\n')
        i = 0
        for ifilename in self.filenames:
            ofile.write("%s %10.8f %12.10f\n" % ( ifilename, self.moduli[i], self.std_moduli[i]) )
            i += 1

        ofile.close()

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument('-f', "--filename", type=str, 
                   help='file with density data.')
    parser.add_argument("-o", "--outfilename", type=str, default='arresttime.txt',
                   help='output file name, assumed to be arresttime.txt')
    parser.add_argument('-l', "--start_line", type=int,
                   help='Set the starting line number.')
        cout << double(i+1) / dN << ' ' << big << ' ' << big*big << ' ' << nclus << ' ' << avg_clus << ' ' << S / dN << '\n';
    parser.add_argument("--occupationfraction_collumn", type=int, default=0,
                   help='Collumn with particle conducting fraction')
    parser.add_argument("--biggestcluster_collumn", type=int, default=1,
                   help='Collumn with the largest cluster size')
    parser.add_argument("--biggestcluster2_collumn", type=int, default=2,
                   help='Collumn with the largest cluster size squared')
    parser.add_argument("--moment_collumn", type=int, default=5,
                   help='Collumn with the second moment')

    # Initialize
    NZ = NewmanZiff()

    err = NZ.addParser(parser)

    if err != None:
        print 'error parsing info'
        return

    NZ.read()
    #G.calc()

    NZ.write()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
