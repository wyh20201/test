#include "matrix.h"

double MAT::GetElem(int h, int l) const
{
    if(h>row-1 || h<0 )
    {
        cout<<"row "<<h<<":"<<"getelem row beyongd range"<<endl;
        exit(1);
    }
    if(l>rank-1 || l<0 )
    {
        cout<<"col"<<l<<":"<<"getelem col beyongd range"<<endl;
        exit(1);
    }
    return this->elem[h+l+h*(this->rank-1)];
}

void MAT::SetElem(int h,int l,double m)
{
    if(h>row-1 || h<0 )
    {
        cout<<"setelem row beyongd range"<<endl;
        exit(1);
    }
    if(l>rank-1 || l<0 )
    {
        cout<<"setelem col beyongd range"<<endl;
        exit(1);
    }
    this->elem[h+l+h*(this->rank-1)]=m;
}

MAT & MAT::operator = (const MAT &other)
{

    if (this == &other)
    {
        return *this;
    }

    delete [] elem; elem=NULL;
    row = other.row;
    rank= other.rank;

    elem = new double[row*rank];
    for(int i=0;i<row*rank;i++)
    {
        elem[i]=other.elem[i];
    }

    return *this;
}

void MAT::Set1()
{
    for (int i = 0; i < row*rank; i++)
    {
        this->elem[i] = 1.0;
    }
}

void MAT::Set0()
{
    for (int i = 0; i < row*rank; i++)
    {
        this->elem[i] = 0.0;
    }
}

void MAT::SetI()
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < rank; j++)
        {
            if (i == j)
            {
                this->SetElem(i, j, 1.0);
            }
            else this->SetElem(i, j, 0.0);
        }
    }
}

MAT::MAT(const MAT &S)
{
    row = S.row;
    rank = S.rank;
    elem = new double[row*rank];
    for (int i = 0; i < row*rank; i++)
    {
        elem[i] = S.elem[i];
    }
}

MAT::MAT(int hang,int lie)
{
    row=hang;
    rank=lie;

    elem=new double[hang*lie];
    for(int i=0;i<row*rank;i++)
    {
        elem[i]=0.0;
    }
}

MAT operator + (const MAT &A,const MAT &B)
{
    int row=A.GetRow();
    int rank=A.GetRank();
    MAT C(row,rank);
    if(row==B.GetRow() && rank==B.GetRank())
    {
        for(int i=0;i<row*rank;i++)
        {
            C.elem[i]=A.elem[i]+B.elem[i];
        }
        return C;
    }
    else
    {
        cout<<"plus matrix error"<<endl;
        exit(1);
    }
}

MAT operator * (const MAT &A,const MAT &B)
{
    int i;
    if((A.GetRank())!=(B.GetRow()))
    {
      cout<<"* matrix error"<<endl;
      exit(1);
    }
    int Crow=A.GetRow();
    int Crank=B.GetRank();
    int Brow=B.GetRow();
    MAT C(Crow,Crank);
    for(i=0;i<(Crow*Crank);i++)
    {
        C.elem[i]=0.0;
    }
    for (i = 0; i < Crow; i++)
    {
        for (int j = 0; j < Crank; j++)
        {
            for (int k = 0; k<Brow; k++)
            {
                C.elem[i + j + i*(Crank - 1)] += A.GetElem(i, k)*B.GetElem(k, j);
            }
        }
    }
        return C;
}

MAT operator * (double A, MAT &B)
{
    int i;
    int Crow = B.GetRow();
    int Crank = B.GetRank();
    MAT C(Crow, Crank);
    for (i = 0; i<(Crow*Crank); i++)
    {
        C.elem[i] = A*B.elem[i];
    }
    return C;
}


void stack::push(int t1,int t2)
{
    node * n=new node(t1,t2,top);
    top=n;
}

node* stack::pop()
{
    node * t=top;
    if(top)
    {
      top=top->prev;
      return t;
    }
    return 0;
}


MAT T(const MAT &B)
{
    int row=B.GetRank();
    int rank=B.GetRow();
    MAT C(row,rank);
    for (int i = 0; i < rank; i++)
    {
        for (int j = 0; j < row; j++)
        {
            int k = i + j + j*(rank - 1);
            C.elem[k] = B.GetElem(i, j);
        }
    }
    return C;
}

int MAT::exrow(int row1,int row2)
{
    if(row1<0 || row1>this->GetRow())
    {
        cout<<"input beyond range"<<endl;
        return 0;
    }
    if(row2<0 || row2>this->GetRow())
    {
        cout<<"input beyond range"<<endl;
        return 0;
    }
    int ran=this->GetRank();
    double* ex;
    ex=new double [ran];
    for (int i = 0; i < ran; i++)
    {
        *(ex + i) = this->GetElem(row1, i);
    }
    for(int i=0;i<ran;i++)
    {
        this->SetElem(row1,i,this->GetElem(row2,i));
        this->SetElem(row2,i,*(ex+i));
    }
    delete [] ex;ex=NULL;
    return 1;
}

int MAT::exrank(int rank1,int rank2)
{
    if(rank1<0 || rank1>this->GetRank())
    {
        cout<<"rank1 beyond range"<<endl;
        return 0;
    }
    if(rank2<0 || rank2>this->GetRank())
    {
        cout<<"rank1 beyond range"<<endl;
        return 0;
    }
    if (rank1 == rank2)
    {
        exit(1);
    }
    int ran=this->GetRow();
    double* ex;
    ex=new double [ran];
    for (int i = 0; i < ran; i++)
    {
        *(ex + i) = this->GetElem(i, rank1);
    }
    for(int i=0;i<ran;i++)
    {
        this->SetElem(i,rank1,this->GetElem(i,rank2));
        this->SetElem(i,rank2,*(ex+i));
    }
    delete [] ex;ex=NULL;
    return 1;
}

MAT inverse1(const MAT &D)
{
    int hang, lie;
    hang = D.GetRow();
    lie = D.GetRank();
    if (hang != lie)
    {
        cout << "row != col";
        exit(2);
    }
    MAT C = D;
    stack exrow;
    for (int i = 0; i<hang; i++)
    {
        int h = i;
        double b = fabs(C.GetElem(i, i));
        if (i<hang - 1)
        {
            for (int j = i; j<hang; j++)
            {
                if (fabs(C.GetElem(j, i))>b)
                {
                    h = j;
                    b = fabs(C.GetElem(j, i));
                }
            }
            if (h != i)
            {
                C.exrow(h, i);
                exrow.push(h, i);
            }
        }
        if (fabs(C.GetElem(i, i))<1e-50)
        {
            cout << "not full";
            exit(0);
        }
        double a = 1 / C.GetElem(i, i);
        C.SetElem(i, i, a);
        for (int j = 0; j < hang; j++)
        {
            if (j != i)
            {
                C.SetElem(j, i, -a*C.GetElem(j, i));
            }
        }
        for (int j = 0; j < hang; j++)
        {
            for (int k = 0; k < hang; k++)
            {
                if (i != k && j != i)
                {
                    C.SetElem(j, k, C.GetElem(j, k) + C.GetElem(j, i)*C.GetElem(i, k));
                }
            }
        }
        for (int j = 0; j < hang; j++)
        {
            if (j != i)
            {
                C.SetElem(i, j, a*C.GetElem(i, j));
            }
        }
    }
    node * p = exrow.top;
    if (p)
    do
    {
        p = exrow.pop();
        C.exrank(p->dat1, p->dat2);
    } while (p->prev != 0);
    return C;
}

bool inverse2(double a[], int ndim, int n, double &d)
{
    int i=0;
    int j=0;
    int k=0;
    int *l; l = new int [ndim];
    int *m; m = new int [ndim];
    double max_a = 0.0;
    double swap  = 0.0;
    d=1.0;
    for(k=0; k<n; k++)
    {
        l[k] = k;
        m[k] = k;
        max_a=a[k*ndim + k];
        for(i=k; i<n; i++)
        {
            for(j=k; j<n; j++)
            {
                if(fabs(max_a)-fabs(a[i*ndim+j]) < 0.0)
                {
                    max_a=a[i*ndim+j];
                    m[k]=i;
                    l[k]=j;
                }
            }
        }
        if(fabs(max_a) < 1.0E-13)
        {
            d=0.0;
            return false;
        }
        if(l[k] > k)
        {
            for(i=0; i<n; i++)
            {
              swap = -a[i*ndim+k];
              a[i*ndim+k]  = a[i*ndim+l[k]];
              a[i*ndim+l[k]] = swap;
            }
        }
        if(m[k] > k)
        {
            for(j=0; j<n; j++)
            {
              swap           = -a[k*ndim+j];
              a[k*ndim+j]    =  a[m[k]*ndim+j];
              a[m[k]*ndim+j] =  swap;
            }
        }
        for(i=0; i<n; i++)
        {
            if(abs(i-k) > 0)
            {
                a[k*ndim+i] = -a[k*ndim+i]/max_a;
            }
        }
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                if((abs(i-k) > 0) && (abs(j-k) > 0))
                {
                    a[j*ndim+i]=a[j*ndim+i]+a[k*ndim+i]*a[j*ndim+k];
                }
            }
        }
        for(j=0; j<n; j++)
        {
            if(abs(j-k) > 0)
            {
                a[j*ndim+k]=a[j*ndim+k]*1.0/max_a;
            }
        }
        a[k*ndim+k]=1.0/max_a;
    }
    for(k=n-1; k>-1; k--)
    {
        if(l[k] > k)
        {
            for(j=0; j<n; j++)
            {
              swap=a[k*ndim+j];
              a[k*ndim+j]=-a[l[k]*ndim+j];
              a[l[k]*ndim+j]=swap;
            }
        }
        if(m[k] > k)
        {
            for(i=0; i<n; i++)
            {
              swap=a[i*ndim+k];
              a[i*ndim+k]=-a[i*ndim+m[k]];
              a[i*ndim+m[k]]=swap;
            }
        }
    }
    delete []l; l=NULL;
    delete []m; m=NULL;
    return true;
}

MAT inverse1new(const MAT &D)
{
    int ndim=D.GetRank();
    double *M2; M2=new double [ndim*ndim];
    for (int i=0;i<ndim;i++)
    {
        for (int j=0;j<ndim;j++)
        {
            M2[i*ndim+j]=D.GetElem(i,j);
        }
    }
    double d=0;
    MAT inv_M2(ndim,ndim);
    inverse2(M2, ndim, ndim,d);
    for (int i=0;i<ndim;i++)
    {
        for (int j=0;j<ndim;j++)
        {
            inv_M2.SetElem(i,j,M2[i*ndim+j]);
        }
    }
    delete []M2; M2=NULL;
    return inv_M2;
}
