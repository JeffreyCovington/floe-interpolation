function  [ee_k,e,fe]=KEspect_2D(u)

Fs=1; %e.g. Fs=1/0.7km samples per km.
 N=min(size(u,1),size(u,2));
 U=u(1:N,1:N);
 dF = Fs/N; f = (-Fs/2:dF:Fs/2-dF); 
   
 func=sin(pi*(0:N-1)/N)'*sin(pi*(0:N-1)/N);

 %e=abs(fftshift(fft2(U.*func))).^2;
 
 e=abs(fftshift(fft2(U))).^2;
 
 [kx, ky]=meshgrid(f,f);
      
   th_int=0:0.001:pi/2; 
   fe=f(f>0);   
  
   [th_int, rho_int]=meshgrid(th_int,fe);
   
   [kx_i, ky_i]=pol2cart(th_int,rho_int);
   
   ee=interp2(kx,ky,e,kx_i, ky_i);
      
   ee_k=fe.*mean(ee');

end