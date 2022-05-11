#include "XNumeric.hpp"
#include "XRandom.hpp"
#include "XDiffEq.hpp"
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>

int N;
int Na;
int Nb;
int P;
int d = 3;
XNum m = 1;
XNum kB = 1;
XNum T;
XNum beta;
XNum hBar = 1.0;
XNum omgP;
XNum L = 30;
XNum omg0 = 1.0 / hBar;
//XNum lamda = 1.0;
XNum s = 0.5;
XNum g = 0.0;
//XNum g = -2.684;
int cor_num = 150;
XNum incre = 0.5*L / cor_num;
XNum vi = 1.0;
/*XNum Qc = 0.0;
XNum sig = 50.0;*/

XNum minimum_image(XNum d) {
  if (std::abs(d) > L / 2)
    return L - std::abs(d);
  return d;
}

XNum minimum_image2(XNum d) {
  if (std::abs(d) > L / 2) {
    if (d < 0)
      return L - std::abs(d);
    else
      return -(L - std::abs(d));
  }
  return d;
}

XNum distance(int i, int j, const std::vector<XNum> &Rs) {
  XNum a = minimum_image(Rs[d*i]-Rs[d*j]);
  XNum b = minimum_image(Rs[d*i+1]-Rs[d*j+1]);
  XNum c = minimum_image(Rs[d*i+2]-Rs[d*j+2]);
  return a*a+b*b+c*c;
}

std::vector<XVecD> ENkCache;

std::pair<XNum, XNum> APhase(const std::vector<XNum> &Rs) {
  std::vector<XNum> RP, IP;
  int N2, k;
  RP.push_back(1);
  IP.push_back(0);
  for (N2 = 1; N2 <= N; ++N2) {
    XNum tmp = 0;
    for (k = 1; k <= N2; ++k) {
      XNum mult = 1;
      //if (((k-1)%2 == 1))
        //mult = -1;
      if (k-1 == 0)
        mult = 1;
      else
        mult = std::pow(vi, k-1);
      tmp += mult*std::exp(-beta*ENkCache[N2-1][k-1])*RP[N2-k];
    }
    RP.push_back(tmp/N2);
    IP.push_back(0);
  }
  std::pair<XNum, XNum> res(RP[N], IP[N]);
  return res;
}

XNum ENk(int N2, int k, const std::vector<XNum> &Rs) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      int index2 = (l - 1) * P + j;
      if (j == P) {
        if (l == N2)
          index2 = (N2 - k) * P;
        else
          index2 = l * P;
      }
      res += distance(index2, index, Rs);
    }
  return 0.5*m*omgP*omgP*res;
}

std::vector<XNum> VBCache;

XNum xexp(XNum k, XNum E, XNum EE) {
  if (vi == 0)
    return std::exp(-beta*E+EE);
  else
    return std::exp((k-1)*std::log(vi)-beta*E+EE);
}

std::vector<XNum> expVB(int N, const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum sum = 0;
    //XNum tmp = 0.5*(ENkCache[N2-1][N2-1]+res[N2-1]);
    XNum tmp = 10000000;
    if (vi == 0)
      tmp = beta*(ENkCache[N2-1][0]+res[N2-1]);
    else {
      for (k = 1; k <= N2; ++k)
        if (-(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+res[N2-k]) < tmp)
          tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+res[N2-k]);
    }
    for (k = 1; k <= N2; ++k) {
      XNum mult = 1;
      /*if (k-1 == 0)
        mult = 1;
      else
        mult = std::pow(vi, k-1);
      sum += mult*std::exp(-beta*(ENkCache[N2-1][k-1]+res[N2-k]-tmp));*/
      if (vi == 0 && k-1 != 0)
        continue;
      sum += mult*xexp(k,ENkCache[N2-1][k-1]+res[N2-k],tmp);
    }
    //res.push_back(tmp-(1/beta)*std::log(sum / N2));
    res.push_back((tmp-std::log(sum)+std::log(N2))/beta);
  }
  return res;
}

std::vector<XNum> gradient(int N2, int k, int l, int j, const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  res.push_back(0);
  res.push_back(0);
  if (l >= N2 - k + 1 && l <= N2) {
    int index = (l - 1) * P + j - 1;
    int index2 = (l - 1) * P + j;
    if (j == P) {
      if (l == N2)
        index2 = (N2 - k) * P;
      else
        index2 = l * P;
    }
    int index3 = (l - 1) * P + j - 2;
    if (j == 1) {
      if (l == N2 - k + 1)
        index3 = (N2 - 1) * P + P - 1;
      else
        index3 = (l - 2) * P + P - 1;
    }
    XNum a = Rs[d*index] - Rs[d*index2];
    XNum b = Rs[d*index] - Rs[d*index3];
    res[0] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
    a = Rs[d*index+1] - Rs[d*index2+1];
    b = Rs[d*index+1] - Rs[d*index3+1];
    res[1] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
    a = Rs[d*index+2] - Rs[d*index2+2];
    b = Rs[d*index+2] - Rs[d*index3+2];
    res[2] = m*omgP*omgP*(minimum_image2(a)+minimum_image2(b));
  }
  return res;
}

std::vector<XNum> XForceVBCache;

void forceVB(int N, const std::vector<XNum> &Rs) {
  int N2, k, l, j;
  XForceVBCache.clear();
  std::vector<XNum> tmp;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      XVecD res;
      res.push_back(0);
      res.push_back(0);
      res.push_back(0);
      for (N2 = 1; N2 <= N; ++N2) {
        int i;
        XNum sum2 = 0;
        //XNum tmp2 = 0.5*(ENkCache[N2-1][N2-1]+VBCache[N2-1]);
        XNum tmp2 = 10000000;
        if (vi == 0)
          tmp2 = beta*(ENkCache[N2-1][0]+VBCache[N2-1]);
        else {
          for (k = 1; k <= N2; ++k)
            if (-(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]) < tmp2)
              tmp2 = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]);
        }
        for (k = 1; k <= N2; ++k) {
          XNum mult = 1;
          /*if (k-1 == 0)
            mult = 1;
          else
            mult = std::pow(vi, k-1);*/
          if (vi == 0 && k-1 != 0)
            continue;
          //sum2 += mult*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
          sum2 += mult*xexp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp2);
        }
        for (i = 0; i < d; ++i) {
          XNum sum = 0;
          for (k = 1; k <= N2; ++k) {
            XNum mult = 1;
            /*if (k-1 == 0)
              mult = 1;
            else
              mult = std::pow(vi, k-1);*/
            if (vi == 0 && k-1 != 0)
              continue;
            tmp = gradient(N2, k, l, j, Rs);
            //sum += mult*(tmp[i] + res[d*(N2-k)+i])*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));
            sum += mult*(tmp[i] + res[d*(N2-k)+i])*xexp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp2);
          }
          res.push_back(sum / sum2);
        }
      }
      XForceVBCache.push_back(res[d*N]);
      XForceVBCache.push_back(res[d*N+1]);
      XForceVBCache.push_back(res[d*N+2]);
    }
}

XVecD force(XNum t, const XVecD &Rs) {
  XVecD res;
  std::vector<XNum> tmp;
  int l, j;
  int tmpN = N;
  N = Na;
  ENkCache.clear();
  XVecD Rsa;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l-1)*P+j-1;
      Rsa.push_back(Rs[d*index]);
      Rsa.push_back(Rs[d*index+1]);
      Rsa.push_back(Rs[d*index+2]);
    }
  if (ENkCache.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j, Rsa);
  VBCache = expVB(N, Rsa);
  forceVB(N, Rsa);
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      XNum a = Rsa[d*index] - L/2;
      XNum b = Rsa[d*index+1] - L/2;
      XNum c = Rsa[d*index+2] - L/2;
      XNum dis = minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b)+minimum_image(c)*minimum_image(c);
      a = Rsa[d*index] - L/2;
      //XNum f = -(1-std::exp(-sig*dis))*Qc*minimum_image2(a)/std::pow(dis,1.5)+2*std::exp(-sig*dis)*Qc*sig*minimum_image2(a)/std::pow(dis,0.5);
      XNum f = 0;
      res.push_back(-XForceVBCache[d*index]/m+f/P-minimum_image2(a)*omg0*omg0/P);
      //res.push_back(-XForceVBCache[d*index]/m-Qc*minimum_image2(a)/(P*std::pow(dis,1.5))+2.0*sig*minimum_image2(a)/(P*std::pow(dis,2.0)));
      a = Rsa[d*index+1] - L/2;
      //f = -(1-std::exp(-sig*dis))*Qc*minimum_image2(a)/std::pow(dis,1.5)+2*std::exp(-sig*dis)*Qc*sig*minimum_image2(a)/std::pow(dis,0.5);
      res.push_back(-XForceVBCache[d*index+1]/m+f/P-minimum_image2(a)*omg0*omg0/P);
      //res.push_back(-XForceVBCache[d*index+1]/m-Qc*minimum_image2(a)/(P*std::pow(dis,1.5))+2.0*sig*minimum_image2(a)/(P*std::pow(dis,2.0)));
      a = Rsa[d*index+2] - L/2;
      //f = -(1-std::exp(-sig*dis))*Qc*minimum_image2(a)/std::pow(dis,1.5)+2*std::exp(-sig*dis)*Qc*sig*minimum_image2(a)/std::pow(dis,0.5);
      res.push_back(-XForceVBCache[d*index+2]/m+f/P-minimum_image2(a)*omg0*omg0/P);
      //res.push_back(-XForceVBCache[d*index+2]/m-Qc*minimum_image2(a)/(P*std::pow(dis,1.5))+2.0*sig*minimum_image2(a)/(P*std::pow(dis,2.0)));
    }
  N = Nb;
  ENkCache.clear();
  Rsa.clear();
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (Na+l-1)*P+j-1;
      Rsa.push_back(Rs[d*index]);
      Rsa.push_back(Rs[d*index+1]);
      Rsa.push_back(Rs[d*index+2]);
    }
  if (ENkCache.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j, Rsa);
  VBCache = expVB(N, Rsa);
  forceVB(N, Rsa);
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l - 1) * P + j - 1;
      XNum a = Rsa[d*index] - L/2;
      XNum b = Rsa[d*index+1] - L/2;
      XNum c = Rsa[d*index+2] - L/2;
      XNum dis = minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b)+minimum_image(c)*minimum_image(c);
      a = Rsa[d*index] - L/2;
      //XNum f = -(1-std::exp(-sig*dis))*Qc*minimum_image2(a)/std::pow(dis,1.5)+2*std::exp(-sig*dis)*Qc*sig*minimum_image2(a)/std::pow(dis,0.5);
      XNum f = 0;
      res.push_back(-XForceVBCache[d*index]/m+f/P-minimum_image2(a)*omg0*omg0/P);
      //res.push_back(-XForceVBCache[d*index]/m-Qc*minimum_image2(a)/(P*std::pow(dis,1.5))+2.0*sig*minimum_image2(a)/(P*std::pow(dis,2.0)));
      a = Rsa[d*index+1] - L/2;
      //f = -(1-std::exp(-sig*dis))*Qc*minimum_image2(a)/std::pow(dis,1.5)+2*std::exp(-sig*dis)*Qc*sig*minimum_image2(a)/std::pow(dis,0.5);
      res.push_back(-XForceVBCache[d*index+1]/m+f/P-minimum_image2(a)*omg0*omg0/P);
      //res.push_back(-XForceVBCache[d*index+1]/m-Qc*minimum_image2(a)/(P*std::pow(dis,1.5))+2.0*sig*minimum_image2(a)/(P*std::pow(dis,2.0)));
      a = Rsa[d*index+2] - L/2;
      //f = -(1-std::exp(-sig*dis))*Qc*minimum_image2(a)/std::pow(dis,1.5)+2*std::exp(-sig*dis)*Qc*sig*minimum_image2(a)/std::pow(dis,0.5);
      res.push_back(-XForceVBCache[d*index+2]/m+f/P-minimum_image2(a)*omg0*omg0/P);
      //res.push_back(-XForceVBCache[d*index+2]/m-Qc*minimum_image2(a)/(P*std::pow(dis,1.5))+2.0*sig*minimum_image2(a)/(P*std::pow(dis,2.0)));
    }
  N = tmpN;
  int k;
  for (j = 1; j <= P; ++j)
    for (l = 1; l <= N; ++l) {
      int index = (l - 1) * P + j - 1;
      int low = 1, high = N;
      for (k = low; k <= high; ++k) {
        if (l == k) continue;
        //if (lamda == 0.0) continue;
        int index2 = (k - 1) * P + j - 1;
        /*XNum inter = lamda / std::pow(distance(index, index2, Rs), 1.5);
        XNum a = Rs[d*index]-Rs[d*index2];
        res[d*index] += minimum_image2(a)*inter/P;
        a = Rs[d*index+1]-Rs[d*index2+1];
        res[d*index+1] += minimum_image2(a)*inter/P;
        a = Rs[d*index+2]-Rs[d*index2+2];
        res[d*index+2] += minimum_image2(a)*inter/P;*/
        XNum inter = (g/(M_PI*s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        XNum a = Rs[d*index]-Rs[d*index2];
        res[d*index] += (2*minimum_image2(a)/(s*s))*inter/P;
        a = Rs[d*index+1]-Rs[d*index2+1];
        res[d*index+1] += (2*minimum_image2(a)/(s*s))*inter/P;
        a = Rs[d*index+2]-Rs[d*index2+2];
        res[d*index+2] += (2*minimum_image2(a)/(s*s))*inter/P;
        /*XNum inter = -(g/(s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        XNum a = Rs[d*index]-Rs[d*index2];
        res[d*index] += (2*minimum_image2(a)/(s*s))*inter/P;
        a = Rs[d*index+1]-Rs[d*index2+1];
        res[d*index+1] += (2*minimum_image2(a)/(s*s))*inter/P;
        a = Rs[d*index+2]-Rs[d*index2+2];
        res[d*index+2] += (2*minimum_image2(a)/(s*s))*inter/P;*/
      }
    }
  return res;
}

XNum GauEnergy(const std::vector<XNum> &Rs) {
  XNum res = 0;
  int j, k, l;
  for (j = 1; j <= P; ++j)
    for (l = 1; l <= N; ++l) {
      int index = (l - 1) * P + j - 1;
      int low = 1, high = N;
      for (k = low; k <= high; ++k) {
        if (l == k) continue;
        //if (lamda == 0.0) continue;
        int index2 = (k - 1) * P + j - 1;
        //XNum inter = lamda / std::pow(distance(index, index2, Rs), 0.5);
        XNum inter = (g/(M_PI*s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        //XNum inter = -(g/(s*s))*std::exp(-distance(index, index2, Rs)/(s*s));
        res += 0.5*inter;
      }
    }
  return res/P;
}

XNum energy(const std::vector<XNum> &Rs) {
  std::vector<XNum> res;
  res.push_back(0);
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    //XNum tmp2 = 0.5*(ENkCache[N2-1][N2-1]+VBCache[N2-1]);
    XNum tmp2 = 10000000;
    if (vi == 0)
      tmp2 = beta*(ENkCache[N2-1][0]+VBCache[N2-1]);
    else {
      for (k = 1; k <= N2; ++k)
        if (-(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]) < tmp2)
          tmp2 = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]);
    }
    XNum sum2 = 0;
    for (k = 1; k <= N2; ++k) {
      XNum mult = 1;
      /*if (k-1 == 0)
        mult = 1;
      else
        mult = std::pow(vi, k-1);
      sum2 += mult*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));*/
      if (vi == 0 && k-1 != 0)
        continue;
      sum2 += mult*xexp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp2);
    }
    XNum sum = 0;
    for (k = 1; k <= N2; ++k) {
      XNum mult = 1;
      /*if (k-1 == 0)
        mult = 1;
      else
        mult = std::pow(vi, k-1);
      sum += mult*(res[N2-k]-ENkCache[N2-1][k-1])*exp(-beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]-tmp2));*/
      if (vi == 0 && k-1 != 0)
        continue;
      sum += mult*(res[N2-k]-ENkCache[N2-1][k-1])*xexp(k,ENkCache[N2-1][k-1]+VBCache[N2-k],tmp2);
    }
    res.push_back(sum / sum2);
  }
  return res[N];
}

XNum energy2(const std::vector<XNum> &Rs) {
  std::vector<XNum> tmp;
  int l, j;
  int tmpN = N;
  N = Na;
  ENkCache.clear();
  XVecD Rsa;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (l-1)*P+j-1;
      Rsa.push_back(Rs[d*index]);
      Rsa.push_back(Rs[d*index+1]);
      Rsa.push_back(Rs[d*index+2]);
    }
  if (ENkCache.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j, Rsa);
  VBCache = expVB(N, Rsa);
  XNum res = energy(Rsa);
  N = Nb;
  ENkCache.clear();
  Rsa.clear();
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = (Na+l-1)*P+j-1;
      Rsa.push_back(Rs[d*index]);
      Rsa.push_back(Rs[d*index+1]);
      Rsa.push_back(Rs[d*index+2]);
    }
  if (ENkCache.empty()) {
    for (l = 1; l <= N; ++l) {
      XVecD tmp;
      for (j = 1; j <= l; ++j)
        tmp.push_back(0);
      ENkCache.push_back(tmp);
    }
  }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j, Rsa);
  VBCache = expVB(N, Rsa);
  res += energy(Rsa);
  N = tmpN;
  return res;
}

void period_boundary(std::vector<XNum> &Rs) {
  int i, j, k;
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      int index = i * P + j;
      for (k = 0; k < d; ++k) {
        if (Rs[d*index+k] < 0) {
          int n = (int)(std::abs(Rs[d*index+k]) / L);
          Rs[d*index+k] += (n + 1) * L;
        }
        else {
          int n = (int)(std::abs(Rs[d*index+k]) / L);
          Rs[d*index+k] -= n * L;
        }
      }
    }
}

XNum box(XNum a) {
  if (a < 0)
    return a+L;
  else if (a >= L)
    return a - L;
  return a;
}

void initial(std::vector<XNum> &Rs, std::vector<XNum> &Vs) {
  int i, j;
  XNum v1 = 0, v2 = 0, v3 = 0;
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      XNum tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v1 += tmp;
      Vs.push_back(tmp);
      tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v2 += tmp;
      Vs.push_back(tmp);
      tmp = XRandGauss() * std::sqrt(1 / (m * beta));
      v3 += tmp;
      Vs.push_back(tmp);
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
      Rs.push_back(L / 2 + 1.0 * (XRandFloat()-0.5));
    }
  v1 /= N * P;
  v2 /= N * P;
  v3 /= N * P;
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      int index = i * P + j;
      Vs[d*index] -= v1;
      Vs[d*index+1] -= v2;
      Vs[d*index+2] -= v3;
    }
}

void test() {
  std::vector<XNum> Rs, Vs;
  int i, j;
  initial(Rs, Vs);
  for (i = 0; i < N; ++i)
    for (j = 0; j < P; ++j) {
      int index = i * P + j;
      std::cout << Rs[index * d] << " " << Rs[index * d + 1] << std::endl;
    }
  XVecD f = force(0, Rs);
  for (i = 0; i < f.size(); ++i)
    std::cout << f[i] << " ";
  std::cout << std::endl;
}

void logR(const std::vector<XNum> &Rs, std::ostream &out) {
  out << "{";
  for (int i = 0; i < N * P; ++i) {
    out << "{" << Rs[d*i] << "," << Rs[d*i+1] << "," << Rs[d*i+2] << "}";
    if (i != N * P - 1)
      out << ",";
  }
  out << "}";
}

void velocity_rescale(std::vector<XNum> &Vs) {
  XNum sum = 0;
  int i;
  for (i = 0; i < N * P; ++i)
    sum += Vs[d*i]*Vs[d*i]+Vs[d*i+1]*Vs[d*i+1]+Vs[d*i+2]*Vs[d*i+2];
  XNum lam = std::sqrt((N*P)*d*(1/beta)/(m*sum));
  for (i = 0; i < N * P; ++i) {
    Vs[d*i] *= lam;
    Vs[d*i+1] *= lam;
    Vs[d*i+2] *= lam;
  }
}

XNum temperature(const std::vector<XNum> &Vs) {
  XNum sum = 0;
  int i;
  for (i = 0; i < N * P; ++i)
    sum += Vs[d*i]*Vs[d*i]+Vs[d*i+1]*Vs[d*i+1]+Vs[d*i+2]*Vs[d*i+2];
  return m*sum/(d*(N*P));
}

XNum trapEnergy(const std::vector<XNum> &Rs) {
  XNum res = 0;
  int i;
  for (i = 0; i < N * P; ++i) {
    XNum a = Rs[d*i] - L/2;
    XNum b = Rs[d*i+1] - L/2;
    XNum c = Rs[d*i+2] - L/2;
    XNum dis = minimum_image(a)*minimum_image(a)+minimum_image(b)*minimum_image(b)+minimum_image(c)*minimum_image(c);
    //res += -Qc/std::pow(dis,0.5)+sig/std::pow(dis,1.0);
    //res += -Qc*(1-std::exp(-sig*dis))/std::pow(dis,0.5);
    res += 0.5*omg0*omg0*dis;
  }
  return res/P;
}

void XMVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, std::vector<XVecD> &theta, std::vector<XVecD> &vtheta, XNum beta, XNum m, std::vector<XVecD> &Q, XForceFunc *func) {
  int i, j;
  int N2 = N;
  int N = d;
  int M = vtheta[0].size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  for (j = 0; j < N2*P; ++j) {
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f2 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    int M2 = M / 2;
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
    for (i = 0; i < N; ++i)
      x[d*j+i] = x[d*j+i] + h*v[d*j+i];
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-1] = theta[j][2*i-1]+h*vtheta[j][2*i-1];
    for (i = 0; i < d; ++i)
      v2[i] = v[d*j+i];
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-2] = vtheta[j][2*i-2]*std::exp(-h*vtheta[j][2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[j][2*i-1]);
  }
  XForceCache = func(t+h, x);
  for (j = 0; j < N2*P; ++j) {
    int M2 = M / 2;
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
  }
}

void XMMVerletNHChain(XNum t, XNum h, XVecD &x, XVecD &v, std::vector<XVecD> &theta, std::vector<XVecD> &vtheta, XNum beta, XNum m, std::vector<XVecD> &Q, XForceFunc *func) {
  int i, j;
  int N2 = N;
  int d2 = d;
  int d = 1;
  int N = d;
  int M = vtheta[0].size();
  if (XForceCache.empty())
    XForceCache = func(t, x);
  for (j = 0; j < N2*P*d2; ++j) {
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f2 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    int M2 = M / 2;
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f2[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
    for (i = 0; i < N; ++i)
      x[d*j+i] = x[d*j+i] + h*v[d*j+i];
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-1] = theta[j][2*i-1]+h*vtheta[j][2*i-1];
    for (i = 0; i < d; ++i)
      v2[i] = v[d*j+i];
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-2] = vtheta[j][2*i-2]*std::exp(-h*vtheta[j][2*i-1])+h*f3[2*i-2]*std::exp(-0.5*h*vtheta[j][2*i-1]);
  }
  XForceCache = func(t+h, x);
  for (j = 0; j < N2*P*d2; ++j) {
    int M2 = M / 2;
    for (i = 0; i < N; ++i)
      v[d*j+i] = v[d*j+i]*std::exp(-0.5*h*vtheta[j][0])+0.5*h*XForceCache[d*j+i]*std::exp(-0.25*h*vtheta[j][0]);
    for (i = 1; i <= M2; ++i)
      theta[j][2*i-2] = theta[j][2*i-2]+h*vtheta[j][2*i-2]/2;
    XVecD v2;
    for (i = 0; i < d; ++i)
      v2.push_back(v[d*j+i]);
    XVecD f3 = NHForce(beta, m, d, Q[j], v2, vtheta[j]);
    for (i = 1; i <= M2; ++i)
      vtheta[j][2*i-1] = vtheta[j][2*i-1]*std::exp(-0.5*h*((i==M2)?0:vtheta[j][2*i]))+0.5*h*f3[2*i-1]*std::exp(-0.25*h*((i==M2)?0:vtheta[j][2*i]));
  }
}

void pair_correlation3(const std::vector<XNum> &Rs, std::vector<XNum> &corrR, std::vector<XNum> &corrI, XNum phaR, XNum phaI) {
  int n = corrR.size();
  int i, j, k;
  for (k = 0; k < P; ++k)
    for (i = 0; i < N; ++i) {
      int in = i*P+k;
      XNum a = Rs[d*in+2]-L/2;
      XNum dis = std::abs(a);
      int index = int(dis / incre);
      if (index < 0)
        index = 0;
      if (index >= n)
        continue;
      corrR[index] += phaR;
      corrI[index] += phaI;
    }
  /*int n = corrR.size();
  int i, j, k;
  for (k = 0; k < P; ++k)
    for (i = 0; i < N; ++i)
      for (j = i + 1; j < N; ++j) {
        XNum dis = std::sqrt(distance(i*P+k, j*P+k, Rs));
        int index = int(dis / incre);
        if (index < 0)
          index = 0;
        if (index >= n)
          index = n - 1;
        corrR[index] += phaR;
        corrI[index] += phaI;
      }*/
}

std::ofstream out("data.txt");

XNum simulation() {
  XForceCache.clear();
  ENkCache.clear();
  VBCache.clear();
  std::vector<XNum> Rs, Vs;
  std::vector<XNum> corrR;
  int i, j;
  initial(Rs, Vs);
  /*logR(Rs, std::cout);
  std::cout << std::endl;*/
  velocity_rescale(Vs);
  XNum t = 0, h = 0.005;
  int skip = 100000;
  for (i = 0; i < cor_num; ++i)
    corrR.push_back(0);
  std::vector<XNum> corrI = corrR;
  std::vector<XVecD> theta, vtheta, Q;
  for (i = 0; i < d*N*P; ++i) {
    XVecD tmp;
    for (j = 0; j < 4; ++j) {
      tmp.push_back(1);
      tmp.push_back(1);
    }
    theta.push_back(tmp);
    tmp.clear();
    for (j = 0; j < 4; ++j) {
      tmp.push_back(1);
      tmp.push_back(1);
    }
    vtheta.push_back(tmp);
    tmp.clear();
    for (j = 0; j < 4; ++j) {
      tmp.push_back(0.1);
      tmp.push_back(0.1);
    }
    Q.push_back(tmp);
  }
  for (i = 0; i < skip; ++i) {
    XMMVerletNHChain(t, h, Rs, Vs, theta, vtheta, beta, m, Q, force);
    period_boundary(Rs);
    t += h;
  }
  logR(Rs, std::cout);
  logR(Rs, out);
  std::cout << std::endl;
  out << std::endl;
  int count = 0;
  int steps = 5000000;
  XNum e1 = 0, e2 = 0, s1 = 0, s2 = 0;
  XNum temp = 0;
  XVecD es1, es2, ss1, ss2;
  XVecD vis;
  for (i = 0; i <= 20; ++i) {
    vis.push_back(-1+0.1*i);
  }
  for (i = 0; i < vis.size(); ++i) {
    es1.push_back(0);
    es2.push_back(0);
    ss1.push_back(0);
    ss2.push_back(0);
  }
  for (i = 0; i < steps; ++i) {
    XMMVerletNHChain(t, h, Rs, Vs, theta, vtheta, beta, m, Q, force);
    period_boundary(Rs);
    if (i % 1 == 0) {
      temp += temperature(Vs);
      XNum te = P * d * N / (2 * beta);
      te += energy2(Rs);
      te += trapEnergy(Rs);
      te += GauEnergy(Rs);
      e1 += te;
      /*std::pair<XNum, XNum> pha = APhase(Rs);
      pha.first /= std::exp(-beta*(VBCache[N]));
      pha.second /= std::exp(-beta*(VBCache[N]));
      e1 += pha.first*te;
      e2 += pha.second*te;
      s1 += pha.first;
      s2 += pha.second;
      XNum tvi = vi;
      for (j = 0; j < vis.size(); ++j) {
        vi = vis[j];
        std::pair<XNum, XNum> pha = APhase(Rs);
        pha.first /= std::exp(-beta*(VBCache[N]));
        pha.second /= std::exp(-beta*(VBCache[N]));
        es1[j] += pha.first*te;
        es2[j] += pha.second*te;
        ss1[j] += pha.first;
        ss2[j] += pha.second;
      }
      vi = tvi;*/
      pair_correlation3(Rs, corrR, corrI, 1, 0);
      ++count;
    }
    if (i % 100000 == 0) {
      std::cout << "i=" << i << std::endl;
      out << "i=" << i << std::endl;
    }
    t += h;
  }
  temp /= count;
  std::cout << temp << std::endl;
  out << temp << std::endl;
  /*e1 /= count;
  e2 /= count;
  s1 /= count;
  s2 /= count;
  std::cout << e1 << "+" << e2 << "i" << std::endl;
  out << e1 << "+" << e2 << "i" << std::endl;
  std::cout << s1 << "+" << s2 << "i" << std::endl;
  out << s1 << "+" << s2 << "i" << std::endl;
  for (j = 0; j < vis.size(); ++j) {
    es1[j] /= count;
    es2[j] /= count;
    ss1[j] /= count;
    ss2[j] /= count;
    std::cout << "vi=" << vis[j] << std::endl;
    std::cout << es1[j] << "+" << es2[j] << "i" << std::endl;
    std::cout << ss1[j] << "+" << ss2[j] << "i" << std::endl;
    out << "vi=" << vis[j] << std::endl;
    out << es1[j] << "+" << es2[j] << "i" << std::endl;
    out << ss1[j] << "+" << ss2[j] << "i" << std::endl;
  }*/
  XNum norm = 0;
  XVecD denR;
  for (i = 0; i < cor_num; ++i) {
    //XNum r = (corrR[i]*s1+corrI[i]*s2)/(s1*s1+s2*s2);
    XNum r = corrR[i];
    denR.push_back(r);
  }
  norm = 0;
  for (i = 0; i < cor_num; ++i)
    norm += denR[i]*incre;
  for (i = 0; i < cor_num; ++i)
    denR[i] /= norm;
  //for (i = 0; i < cor_num; ++i)
    //denR[i] /= 2*M_PI*(i+0.5)*incre;
  for (i = 0; i < cor_num; ++i)
    std::cout << "{" << i*incre << "," << denR[i] << "},";
  std::cout << std::endl;
  for (i = 0; i < cor_num; ++i)
    out << "{" << i*incre << "," << denR[i] << "},";
  out << std::endl;
  logR(Rs, std::cout);
  logR(Rs, out);
  std::cout << std::endl;
  out << std::endl;
  /*return e1;*/
  return e1/count;
}

int main(int argc, char *argv[]) {
  clock_t t;
  t = clock();
  Na = 5;
  Nb = 5;
  N = Na + Nb;
  T = 1.0 / 1.0;
  beta = 1 / (kB * T);
  /*for (P = 12; P <= 12; P += 3) {
    omgP = std::sqrt(P) / (beta * hBar);
    XNum e = simulation();
    std::cout << "{" << P << ", " << e << "}, ";
    out << "{" << P << ", " << e << "}, ";
  }*/
  XNum vv = 0;
  if (argc >= 2) {
    out.close();
    out.open(std::string(argv[1]));
  }
  if (argc >= 3)
    vv = std::atof(argv[2]);
  for (vi = vv; vi <= vv; vi += 0.1) {
    XSetRandSeed(5);
    P = 50;
    omgP = std::sqrt(P) / (beta * hBar);
    XNum e = simulation();
    std::cout << "{" << vi << ", " << e << "}, ";
    out << "{" << vi << ", " << e << "}, ";
    std::cout << std::endl;
    out << std::endl;
  }
  std::cout << std::endl;
  out << std::endl;
  t = clock() - t;
  std::cout << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  out << (int)(((double)1000 * t) / CLOCKS_PER_SEC) << std::endl;
  return 0;
}
