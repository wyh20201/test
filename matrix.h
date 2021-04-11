#ifndef MATRIX_H
#define MATRIX_H
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
using namespace std;
class MAT
{

public:
    MAT()
    {
        row=rank=0;
        elem=new double;
        *elem=0.0;
    }
    MAT(int hang,int lie);
    MAT(const MAT &S);
    ~MAT(){ delete [] elem;elem=NULL;}
    void SetI();
    void Set1();
    void Set0();
    void SetElem(int h,int l,double m);
    void SetRow(int h)
    {
        row=h;
        SetRR();
    }
    void SetRank(int l)
    {
        rank=l;
        SetRR();
    }
    int GetRow() const
    {
        return this->row;
    }
    int GetRank() const
    {
        return this->rank;
    }
    double GetElem(int h, int l)const;
    int exrow(int,int );
    int exrank(int,int);
    MAT & operator = (const MAT &other);
    friend MAT operator + (const MAT &A,const MAT &B);
    friend MAT operator * (const MAT &A,const MAT &B);
    friend MAT operator * (double A, MAT &B);
    friend MAT T(const MAT &B);
    friend MAT inverse1(const MAT &D);
    friend bool inverse2(double a[], int ndim, int n, double &d);
    friend MAT inverse1new(const MAT &D);

private:
    int row;
    int rank;
    double *elem;
    void SetRR()
    {
        delete [] elem; elem=NULL;
        this->elem = new double[GetRow()*GetRank()];
    }
};

class node
{
   public:
        int dat1;
        int dat2;
        node * prev;
        node(int d1,int d2,node *n)
        {
             dat1=d1;
             dat2=d2;
             prev=n;
        }
        ~node(){}
};

class stack
{
    public:
        node *top;
        stack()
        {
            top=0;
        }
        void push(int,int);
        node* pop();
};
#endif
