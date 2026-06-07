#ifndef __KD_PHOTON_H__
#define __KD_PHOTON_H__
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;

#define KD_HEADER_ID 0xABCD
#define KD_PHOTON_ID 0x0001

typedef  unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned long uint32;



struct KD_HEADER
{                               //81-336B
    uint16 hid;                 //2B
    uint16 hsize;               //2B
    uint8 char_num;             //1B
    char descr[255];            //0-255B
    double a;                   //8B
    double rin, rout;           //16B
    double ro,tho;              //16B
    double alphamin,alphamax;   //16B
    double betamin,betamax;     //16B
    uint16 nalpha,nbeta;        //4B
    uint32 nphotons;            //4B
};

struct KD_PHOTON
{                               //42+8*ni
    uint16 hid;                 //2B
    uint16 hsize;               //2B
    double l,q;                 //16B
    double alpha,beta;          //16B
    uint16 i,j;                    //4B
    uint16 num_normal;             //2B
    uint16 num_reflx;              //2B
    vector<double> re_n;        //num_normal*8B
    vector<double> re_r;        //num_reflx*8B
};

#endif
