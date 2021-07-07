/*
Calculates the largest site and bond percolated cluster using the Newman-Ziff method
Newman, M. E. J., & Ziff, R. M. (2001). Fast Monte Carlo algorithm for site or bond percolation. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, 64(1), 16. https://doi.org/10.1103/PhysRevE.64.016706
*/
#include "newmanziff.h"

int newmanziff::get_NN(const std::string& forfile, const int& max_NN, const int& npoints, int nsteps) {

    std::istringstream ssf;
   
    int tstepF = 0;
    int ncontacts;
  
    // find final timestep
    if (nsteps < 0) {
        std::ifstream tfin(forfile);
        std::string tfline;
        int nsteps = 0;
        while (std::getline(tfin,tfline)) if (tfline=="ITEM: TIMESTEP") nsteps++;
        tfin.close();
    } else if (nsteps == 0) {
        std::cout << "Frames are indexed at 1 not 0\n";
        return -1;
    }

    std::ifstream fin(forfile);
    std::string fline;

    int istep = 1;
    while (std::getline(fin,fline)) {
        if (fline=="ITEM: TIMESTEP") {
            if (istep == nsteps) {
                // read header
                std::getline(fin,fline); 
                ssf.str(fline); 
                ssf >> tstepF;
                ssf.clear();
                std::getline(fin,fline);
                std::getline(fin,fline);
                ssf.str(fline);
                // total number of contacts between partcles or bonds between sites
                ssf >> ncontacts;
                ssf.clear();

                //int nn[ncontacts][max_NN]; /* Nearest neighbors */
                int nc[npoints]; /* Nearest neighbors */
                nn = new int*[npoints];
                for (int i = 0; i < npoints; i++) {
                    nc[i] = 0;
                    nn[i] = new int[max_NN];
                    for (int j = 0; j < max_NN; j++) {
                        nn[i][j] = -1;
                    }
                }
                for (int iline = 0; iline < 5; iline++) { std::getline(fin,fline); }
  
                int id, iatom, jatom;
                for (int i=0;i<ncontacts;i++) {
                    std::getline(fin,fline);
                    ssf.str(fline);
                    ssf >> id >> iatom >> jatom;
                    nn[iatom-1][nc[iatom-1]++] = jatom-1;
                    ssf.clear();
                }
                for (int iatom = 0; iatom < npoints; iatom++) {
                    if ( nc[iatom] > max_NN ) {
                        std::cout << "Atom " << iatom << " has " << nc[iatom] << " contacts, you need to increase max_NN above " << max_NN << "\n";
                        return -1;
                    }
                }
            }
            istep++;
        }
    }

    if (istep < nsteps+1) {
        std::cout << nsteps << " frames requested, but there are only " << istep << " frames\n";
        return -1;
    }
    fin.close();
    return 0;
}

void newmanziff::permutation( const int& N )
{
    order = new int[N]; /* Occupation order */
    int i,j;
    int temp;
    for (i=0; i<N; i++) order[i] = i;
    for (i=0; i<N; i++) {
        j = i + (N-i) * double(rand()) / RAND_MAX;
        temp = order[i];
        order[i] = order[j];
        order[j] = temp;
    }
}

int newmanziff::findroot(int i)
{
    if (ptr[i]<0) return i;
    return ptr[i] = findroot(ptr[i]);
}

void newmanziff::percolate( const std::string& outfile, const int& N, const int& max_NN )
{
    int EMPTY = -N-1; 

    std::ostringstream oss;
    oss << "# fraction_occupied max_cluster max_cluster^2 avg_cluster/N S\n";
    std::string header = oss.str();
    std::ofstream cout(outfile);
    cout << header;

    ptr = new int[N]; /* Array of pointers */

    int i,j;
    int s1,s2; // site index
    int r1,r2;
    int big=0;
    double dN = double(N);
    int nclus=0;
    double avg_clus;
    int S = 0; // second moment
    for (i=0; i<N; i++) ptr[i] = EMPTY;
    for (i=0; i<N; i++) {
        r1 = s1 = order[i];
        ptr[s1] = -1;
        int bound = 0;
        nclus++;
        for (j=0; j<max_NN; j++) {
            s2 = nn[s1][j]; // site index for this neighbor j

            if (s2 < 0) break; // reached the limit of contacts for this particle

            // neighbor is already part of a cluster
            if (ptr[s2]!=EMPTY) {
                bound = 1;
                r2 = findroot(s2);  // find the root cluster
                if (r2!=r1) {  // are they as of now, in different clusters
                    nclus--;  // combine clusters
                    S += 2 * ptr[r1] * ptr[r2]; // add 2 * 
                    if (ptr[r1]>ptr[r2]) {
                        ptr[r2] += ptr[r1];
                        ptr[r1] = r2;
                        r1 = r2;
                    } else {
                        ptr[r1] += ptr[r2];
                        ptr[r2] = r1;
                        r2 = r1;
                    }
                    if (-ptr[r1]>big) big = -ptr[r1];
                }
            }
        }
        // add a cluster of size 1 since it was not connected
        if ( bound == 0 ) {
            S++; // 1^2= 1
        }
        if (nclus < 1) {
            avg_clus = i+1;
        } else {
            avg_clus = double(i+1)/double(nclus);
        }
        cout << double(i+1) / dN << ' ' << big << ' ' << big*big << ' ' << nclus << ' ' << avg_clus << ' ' << S / dN << '\n';
    }
  cout.close();
}


