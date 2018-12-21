
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include "../include/fileIO.h"

namespace FileIO{

    /*
     * Reads datafile into memory.
     */
    double2* readIn(std::string fileR, std::string fileI,
                        int gSize){
        FILE *f;
        f = fopen(fileR.c_str(),"r");
        int i = 0;
        double2 *arr = (double2*) malloc(sizeof(double2)*gSize);
        double line;
        while(fscanf(f,"%lE",&line) > 0){
            arr[i].x = line;
            ++i;
        }
        fclose(f);
        f = fopen(fileI.c_str(),"r");
        i = 0;
        while(fscanf(f,"%lE",&line) > 0){
            arr[i].y = line;
            ++i;
        }
        fclose(f);
        return arr;
    }

    /*
     * Writes out the parameter file.
     */
    void writeOutParam(Grid &par, std::string file){
        par.write(file);
    }

    /*
     * Writes out double2 complex data files.
     */
    void writeOut(std::string buffer, std::string file, double2 *data,
                      int length, int step){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_%d", file.c_str(), step);
        f = fopen (buffer.c_str(),"w");
        int i;
        for (i = 0; i < length; i++)
            fprintf (f, "%.16e\n",data[i].x);
        fclose (f);

        sprintf ((char *)buffer.c_str(), "%si_%d", file.c_str(), step);
        f = fopen (buffer.c_str(),"w");
        for (i = 0; i < length; i++)
            fprintf (f, "%.16e\n",data[i].y);
        fclose (f);

    }

    /*
     * Writes out double type data files.
     */
    void writeOutDouble(std::string file, double *data, int length, int step){
        std::ofstream output;
        output.open(file + "_" + std::to_string(step));
        for (int i = 0; i < length; ++i){
            output << data[i] << '\n';
        }

        output.close();
    }

    /*
     * Writes out bool type data files.
     */
    void writeOutBool(std::string file, bool *data,int length, int step){
        std::ofstream output;
        output.open(file + "_" + std::to_string(step));
        for (int i = 0; i < length; ++i){
            output << data[i] << '\n';
        }

        output.close();
    }

    /*
     * Writes out int type data files.
     */
    void writeOutInt(std::string buffer, std::string file, int *data,
                         int length, int step){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_%d", file.c_str(), step);
        f = fopen (buffer.c_str(),"w");
        int i;
        for (i = 0; i < length; i++)
            fprintf (f, "%d\n",data[i]);
        fclose (f);
    }

    /*
     * Writes out int2 data type.
     */
    void writeOutInt2(std::string buffer, std::string file, int2 *data,
                          int length, int step){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_%d", file.c_str(), step);
        f = fopen (buffer.c_str(),"w");
        int i;
        for (i = 0; i < length; i++)
            fprintf (f, "%d,%d\n",data[i].x,data[i].y);
        fclose (f);
    }

    /*
     * Writes out tracked vortex data.
     */
    void writeOutVortex(std::string buffer, std::string file,
                            std::vector<std::shared_ptr<Vtx::Vortex>> &data, int step){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_%d", file.c_str(), step);

        f = fopen (buffer.c_str(),"w");
        int i;

        fprintf (f, "#UID,X,Xd,Y,Yd,WINDING,isOn\n");
        for (i = 0; i < data.size(); i++)
            //fprintf (f, "%d,%d,%e,%d,%e,%d\n",data[i]->getUID(),data[i]->getCoords().x,data[i]->getCoordsD().x,data[i]->getCoords().y,data[i]->getCoordsD().y,data[i]->getWinding());
            fprintf (f, "%d,%e,%d,%e,%d\n",data[i]->getCoords().x,data[i]->getCoordsD().x,data[i]->getCoords().y,data[i]->getCoordsD().y,data[i]->getWinding());
        fclose (f);
    }

    /*
     * Opens and closes file. Nothing more. Nothing less.
     */
    int readState(std::string name){
        FILE *f;
        f = fopen(name.c_str(),"r");
        fclose(f);
        return 0;
    }

    /*
     * Outputs the adjacency matrix to a file
     */
    void writeOutAdjMat(std::string buffer, std::string file, int *mat, unsigned int *uids, int dim, int step){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_%d", file.c_str(), step);
        f = fopen (buffer.c_str(),"w");
        fprintf (f, "(*");
        for(int ii = 0; ii<dim; ++ii){
            fprintf (f, "%d",uids[ii]);
        }
        fprintf (f, "*)\n");
        fprintf (f, "{\n");
        for(int ii = 0; ii < dim; ++ii){
            fprintf (f, "{");
            for(int jj = 0; jj < dim; ++jj){
                fprintf (f, "%i",mat[ii*dim + jj]);
                if(jj<dim-1)
                    fprintf (f, ",");
                else
                    fprintf (f, "}");
            }
            if(ii<dim-1)
                fprintf (f, ",");
            fprintf (f, "\n");
        }
        fprintf (f, "}\n");
        fclose(f);
    }
    void writeOutAdjMat(std::string buffer, std::string file, double *mat,
                        unsigned int *uids, int dim, int step){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_%d", file.c_str(), step);
        f = fopen (buffer.c_str(),"w");
        fprintf (f, "(*");
        for(int ii = 0; ii<dim; ++ii){
            fprintf (f, "%d",uids[ii]);
            if(ii!=dim-1)
               /* I am not sure what Lee wants here, but I think...
                           fprintf (f, ",",uids[ii]); */
                           fprintf (f, ",");

        }
        fprintf (f, "*)\n");
        fprintf (f, "{\n");
        for(int ii = 0; ii < dim; ++ii){
            fprintf (f, "{");
            for(int jj = 0; jj < dim; ++jj){
                fprintf (f, "%e",mat[ii*dim + jj]);
                if(jj<dim-1)
                    fprintf (f, ",");
                else
                    fprintf (f, "}");
            }
            if(ii<dim-1)
                fprintf (f, ",");
            fprintf (f, "\n");
        }
        fprintf (f, "}\n");
        fclose(f);
    }
}
