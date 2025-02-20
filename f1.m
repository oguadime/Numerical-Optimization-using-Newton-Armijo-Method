function [f,fp,fpp]=f1(x) 
  global alpha;
  
  f = (sin(x))^2+alpha*rand(1);
  fp = 2*sin(x)*cos(x);
  fpp = 2*((cos(x))^2-(sin(x))^2);
