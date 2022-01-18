/*
Calculates the largest site and bond percolated cluster using the Newman-Ziff method
Newman, M. E. J., & Ziff, R. M. (2001). Fast Monte Carlo algorithm for site or bond percolation. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics, 64(1), 16. https://doi.org/10.1103/PhysRevE.64.016706
*/
#include "newmanziff.h"

void newmanziff::remove_rattlers(std::vector<int>& i_contacts, std::vector<int>& j_contacts, std::vector<int>& ncont_atom, const float& local_isostatic) {

  int ncont = i_contacts.size();
    
  for (int i=0;i<ncont;i++) {
    ncont_atom[i_contacts[i]-1] += 1;
    ncont_atom[j_contacts[i]-1] += 1;
  }
  
  std::vector<int> i_contacts_new; i_contacts_new.reserve(i_contacts.size());
  std::vector<int> j_contacts_new; j_contacts_new.reserve(j_contacts.size());
  for (int i=0;i<ncont;i++) {
    if ((ncont_atom[i_contacts[i]-1] > local_isostatic) && (ncont_atom[j_contacts[i]-1] > local_isostatic)) {
      i_contacts_new.emplace_back(i_contacts[i]);    
      j_contacts_new.emplace_back(j_contacts[i]);
    }
  }

  int ncont_new = i_contacts_new.size();
  i_contacts.assign(ncont_new,0); 
  i_contacts.shrink_to_fit();
  j_contacts.shrink_to_fit();
  for (int i=0;i<ncont_new;i++) {
    i_contacts[i] = i_contacts_new[i];
    j_contacts[i] = j_contacts_new[i];
  }
}
   
int newmanziff::get_NN(const std::string& forfile, const int& max_NN, const int& npoints, int nsteps, const int& Nrealizations, const float& local_isostatic) {

    std::istringstream ssf;
   
    int tstepF = 0;
    int ncontacts;
  
    // find final timestep
    if (nsteps < 0) {
        std::ifstream tfin(forfile);
        std::string tfline;
        int isteps = 0;
        while (std::getline(tfin,tfline)) if (tfline=="ITEM: TIMESTEP") isteps++;
        tfin.close();
        nsteps = isteps;
    } else if (nsteps == 0) {
        std::cout << "Frames are indexed at 1 not 0\n";
        return -1;
    }
    std::cout << nsteps << "nsteps\n";

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
                for (int iline = 0; iline < 5; iline++) { std::getline(fin,fline); }
  
                std::vector<int> contacts_i(ncontacts);
                std::vector<int> contacts_j(ncontacts);
                int id, iatom, jatom;
                for (int i=0;i<ncontacts;i++) {
                    std::getline(fin,fline);
                    ssf.str(fline);
                    ssf >> id >> iatom >> jatom; // >> fn  >> fs >> fr >> ft;
                    contacts_i[i] = iatom;
                    contacts_j[i] = jatom;
                    //nn[iatom-1][nc[iatom-1]++] = jatom-1;
                    ssf.clear();
                }

                std::vector<int> conti_temp(ncontacts);
                std::vector<int> contj_temp(ncontacts);
              
                for (int i=0;i<ncontacts;i++) {
                  conti_temp[i] = contacts_i[i];
                  contj_temp[i] = contacts_j[i];
                }
                
                int nc_new = ncontacts;

                // remove rattlers
                if (local_isostatic > 0) {
                  std::vector<int> ncont_atom(npoints,0);  // number of contacts at each atom
                  while (1) {
                    int tnc = conti_temp.size();
                    remove_rattlers(conti_temp, contj_temp, ncont_atom, local_isostatic);
                    nc_new = conti_temp.size();
                    if ((nc_new == tnc) || (nc_new==0)) {
                      break;
                    }
                  }
                }

                int nc[npoints]; /* Nearest neighbors */
                nn = new int*[npoints];
                for (int i = 0; i < npoints; i++) {
                    nc[i] = 0;
                    nn[i] = new int[max_NN];
                    for (int j = 0; j < max_NN; j++) {
                        nn[i][j] = -1;
                    }
                }
/*
                for (int i=0;i<ncontacts;i++) {
                    nn[contacts_i[i]-1][nc[contacts_i[i]-1]++] = contacts_j[i]-1;
                    ssf.clear();
                }
*/
                for (int i=0;i<nc_new;i++) {
                    nn[conti_temp[i]-1][nc[conti_temp[i]-1]++] = contj_temp[i]-1;
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

    // initialize arrays
    avg_cluster = new double*[Nrealizations]; /* Array of pointers */
    max_cluster = new int*[Nrealizations]; /* Array of pointers */
    S = new int*[Nrealizations]; /* Array of pointers */
    for (int i = 0; i < Nrealizations; i++) {
        avg_cluster[i] = new double[npoints];
        max_cluster[i] = new int[npoints];
        S[i] = new int[npoints];
    }

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

//void newmanziff::percolate( const std::string& outfile, const int& N, const int& max_NN )
void newmanziff::percolate( const int& N, const int& max_NN, const int& realization )
{
    int EMPTY = -N-1; 

/*
    std::ostringstream oss;
    oss << "# fraction_occupied max_cluster max_cluster^2 avg_cluster/N S\n";
    std::string header = oss.str();
    std::ofstream cout(outfile);
    cout << header;
*/

    ptr = new int[N]; /* Array of pointers */

    int i,j;
    int s1,s2; // site index
    int r1,r2;
    int big = 0;
    double dN = double(N);
    int nclus = 0;
    //double avg_clus;
    int iS = 0; // second moment
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
                    iS += 2 * ptr[r1] * ptr[r2]; // add 2 * 
                    //S[realization][i] += 2 * ptr[r1] * ptr[r2]; // add 2 * 
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
            iS++; // 1^2= 1
        }
        if (nclus < 1) {
            avg_cluster[realization][i] = double(i+1);
        } else {
            avg_cluster[realization][i] = double(i+1)/double(nclus);
        }
        S[realization][i] = iS / dN;
        max_cluster[realization][i] = big;
        //std::cout << double(i+1) / dN << ' ' << max_cluster[realization][i] << ' ' << max_cluster[realization][i]*max_cluster[realization][i] << ' ' << nclus << ' ' << avg_cluster[realization][i] << ' ' << S[realization][i] / dN << '\n';
        //std::cout << double(i+1) / dN << ' ' << big << ' ' << big*big << ' ' << nclus << ' ' << avg_cluster[realization][i] << ' ' << iS / dN << '\n';
    }
  //cout.close();
}


void newmanziff::write( const std::string& outfile, const int& N, const int& Nrealizations )
{
    double dN = double(N);
    double dNr = double(Nrealizations);
    // calculate average and standard deviations
    // calculate average and standard deviations
    double allbig[N];
    double allbigbig[N];
    double allave[N];
    double allS[N];
    double stdbig[N];
    double stdbigbig[N];
    double stdave[N];
    double stdS[N];

    for (int i=0; i<N; i++) {
        allbig[i] = 0.0;
        allbigbig[i] = 0.0;
        allave[i] = 0.0;
        allS[i] = 0.0;
        stdbig[i] = 0.0;
        stdbigbig[i] = 0.0;
        stdave[i] = 0.0;
        stdS[i] = 0.0;
        for (int r=0; r<Nrealizations; r++) {
            //std::cout << "i r " << i << ' ' << r << ' ' << max_cluster[r][i] << ' ' << avg_cluster[r][i] << ' ' << S[r][i] <<'\n';
            allbig[i] += double(max_cluster[r][i]);
            allbigbig[i] += double(max_cluster[r][i]) * double(max_cluster[r][i]);
            allave[i] += avg_cluster[r][i];
            allS[i] += S[r][i];
        }
        allbig[i] /= dNr;
        allbigbig[i] /= dNr;
        allave[i] /= dNr;
        allS[i] /= dNr;
        for (int r=0; r<Nrealizations; r++) {
            double d1 = (max_cluster[r][i] - int(allbig[i]));
            double d2 = ((max_cluster[r][i]*max_cluster[r][i]) - int(allbigbig[i]));
            double d3 = (avg_cluster[r][i] - int(allave[i]));
            double d4 = (S[r][i] - int(allS[i]));
            stdbig[i] += d1 * d1;
            stdbigbig[i] += d2 * d2;
            stdave[i] += d3 * d3;
            stdS[i] += d4 * d4;
        }
        stdbig[i] = sqrt(stdbig[i] / dNr);
        stdbigbig[i] = sqrt(stdbigbig[i] / dNr);
        stdave[i] = sqrt(stdave[i] / dNr);
        stdS[i] = sqrt(stdS[i] / dNr);
    }

    std::ostringstream oss;
    oss << "# fraction_occupied <max_cluster> <max_cluster^2> <avg_cluster> <moment> std(max_cluster) std(max_cluster^2) std(avg_cluster) std(moment)\n";
    std::string header = oss.str();
    std::ofstream cout(outfile);
    cout << header;

    for (int i=0; i<N; i++) {
        cout << double(i+1) / dN << ' ' << allbig[i] << ' ' << allbigbig[i] << ' ' << allave[i] << ' ' << allS[i] << ' ' << stdbig[i] << ' ' << stdbigbig[i] << ' ' << stdave[i] << ' ' << stdS[i] << '\n'; 
    }
  cout.close();
}
