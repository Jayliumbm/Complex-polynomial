#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

class C{//complex number

public:
    double re,im;
    C()
        :re(0), im(0){}

    C(double r, double i)
        :re(r), im(i){}

    friend ostream& operator<<(ostream& os, const C& a);

    friend istream& operator<<(istream& in, C& a);

    double mod(){
        return sqrt(re*re+im*im);
    }
    double mod2(){
        return (re*re+im*im);
    }

    C operator+(C b){
        C c;
        c.re=re+b.re;
        c.im=im+b.im;
        return c;
    }

    C& operator+=(C b){
        this->re+=b.re;
        this->im+=b.im;
        return *this;
    }

    C operator-(C b){
        C c;
        c.re=re-b.re;
        c.im=im-b.im;
        return c;
    }

    C& operator-=(const C b){
        this->re-=b.re;
        this->im-=b.im;
        return *this;
    }

    C operator*(C b){
        C c;
        c.re=re*(b.re)-im*(b.im);
        c.im=im*(b.re)+re*(b.im);
        return c;
    }

    C operator*(double r){
        C c(re*r,im*r);
        return c;
    }

    C& operator*=(const C b){
        double r=this->re;
        double i=this->im;
        this->re=r*(b.re)-i*(b.im);
        this->im=i*(b.re)+r*(b.im);
        return *this;
    }

    C& operator*=(double r){
        this->re*=r;
        this->im*=r;
        return *this;
    }

    C operator/(C b){
        C c;
        c.re=(re*(b.re)+im*(b.im))/b.mod2();
        c.im=(im*(b.re)-re*(b.im))/b.mod2();
        return c;
    }

    C operator/(double r){
        C c(re/r,im/r);
        return c;
    }

    C& operator/=(C b){
        double r=this->re;
        double i=this->im;
        this->re=(r*(b.re)+i*(b.im))/b.mod2();
        this->im=(i*(b.re)-r*(b.im))/b.mod2();
        return *this;
    }

    C& operator/=(double r){
        this->re/=r;
        this->im/=r;
        return *this;
    }

    void set_val(double r,double i){
        re=r;
        im=i;
    }
};

ostream& operator<<(ostream& os, const C& a){

    double r=(abs(a.re)<=0.00049)?0:a.re;
    double i=(abs(a.im)<=0.00049)?0:a.im;

    os<<fixed<<setprecision(3)<<((r>=0)?" + ":" - ")<<abs(r)<<((i>=0)?" + ":" - ")<<abs(i)<<" i";

    return os;
}

istream& operator>>(istream& in, C& a){

    in >> a.re;
    in >> a.im;
    return in;
}

class func{
public:
    int deg;
    vector<C> coif;
    friend ostream& operator<<(ostream& os, const func& a);
    func(vector<C> co)
        :coif(co)
    {
        deg=co.size()-1;
    };//constructor

    C fv(C);//find value
};

ostream& operator<<(ostream& os, const func& a){
    for(int i=0;i<=a.deg;i++){

        os<<a.coif[i]<<endl;
    }
    return os;
}

C func::fv(C z){

    C fz;

    for(int i=0;i<=deg;i++){
        C zn(1,0);
        for(int j=1;j<=i;j++){
            zn*=z;
        }
        fz+=(zn*coif[i]);
    }
    return fz;
}

int root_num(func f, func df, C center, double radius){

    double x1=center.re-radius;
    double x2=center.re+radius;
    double y1=center.im-radius;
    double y2=center.im+radius;
    double dx=radius/1000,dy=radius/1000;

    C I;
    C z(0,0);

    C dz(-dx,0);
    while(x1<x2){//straight line from (x2,y2) to (x1,y2)
        z.set_val(x2,y2);
        I+=(df.fv(z)/f.fv(z)*dz);
        x2-=dx;
    }
    x2=center.re+radius;

    dz.set_val(dx,0);
    while(x1<x2){//straight line from (x1,y1) to (x2,y1)
        z.set_val(x1,y1);
        I+=(df.fv(z)/f.fv(z)*dz);
        x1+=dx;
    }
    x1=center.re-radius;

    dz.set_val(0,-dy);
    while(y1<y2){//straight line from (x1,y2) to (x1,y1)
        z.set_val(x1,y2);
        I+=(df.fv(z)/f.fv(z)*dz);
        y2-=dy;
    }
    y2=center.im+radius;

    dz.set_val(0,dy);
    while(y1<y2){//straight line from (x2,y1) to (x2,y2)
        z.set_val(x2,y1);
        I+=(df.fv(z)/f.fv(z)*dz);
        y1+=dy;
    }
    y1=center.im-radius;
    return ((int)(I.im/3.14159/2+0.5));//Imaginary part divided by 2pi
}


vector<C> solve(int d,func f,func df,double range){

    vector<C> roots;
    double radius=range+0.1;
    for(int i=0;i<d;i++){
        C c(0.1,0.1);
        roots.push_back(c);
    }
    C v[4];
    v[0].set_val(0.5,0.5);
    v[1].set_val(-0.5,0.5);
    v[2].set_val(-0.5,-0.5);
    v[3].set_val(0.5,-0.5);

    int done=0;
    while(!done){
        int n[4];
        vector<C> new_roots;
        for(int i=0;i<d;i++){
            for(int dir=0;dir<4;dir++){
                n[dir]=root_num(f,df,roots[i]+v[dir]*radius,0.5*radius);
                for(int j=0;j<n[dir];j++){
                    new_roots.push_back(roots[i]+v[dir]*radius);
                }
            }
            i+=n[0]+n[1]+n[2]+n[3]-1;
        }
        roots=new_roots;
        if(radius<0.001){
            done=1;
        }
        radius/=2;
    }
    return roots;
}


int main () {

    cout<<"This is a calculator for finding roots of complex-coefficient polynomials."<<endl;
    cout<<"Does your polynomial has complex coefficient? If yes, enter 1; if no, enter 0.\n"<<endl;
    int is_complex;
    cin>>is_complex;

    cout<<"\nPlease enter the degree of your polynomial."<<endl;
    cout<<"(less than 5 would be better)\n"<<endl;
    int deg;
    cin>>deg;

    vector<C> co;
    if(is_complex){
        cout<<"\nEnter your polynomial in terms of coefficients in descending power."<<endl;
        cout<<"For example, z^2+(4+i)z+(5-2i)=0 =>\"1 0 4 1 5 -2\" \n"<<endl;

        for(int i=0;i<=deg;i++){
            C c;
            cin>>c;
            co.insert(co.begin(), c);
        }
    }

    else{
        cout<<"\nEnter your polynomial in terms of coefficients in descending power."<<endl;
        cout<<"For example, z^3-2z^2+4z+5=0 =>\"1 -2 4 5\" \n"<<endl;

        for(int i=0;i<=deg;i++){
            double a;
            cin>>a;
            C c(a,0);
            co.insert(co.begin(), c);
        }
    }

    func f(co);

    double range=0;
    for(int i=0;i<f.deg;i++){//find max{|c_0|,|c_1|, ... ,|c_(n-1)|}
        if(f.coif[i].mod()>range){
            range=f.coif[i].mod();
        }
    }

    range=range*deg/f.coif[deg].mod()+1;

    vector<C> dco;

    for(int i=0;i<deg;i++){
        dco.push_back(co[i+1]*(double)(i+1));
    }
    func df(dco);

    vector<C> roots=solve(deg,f,df,range);

    cout<<endl<<"All roots:"<<endl;

    for(int i=0;i<deg;i++){
        cout<<roots[i]<<endl;
    }

    return 0;
}
