tola=1e-6;
tolr=tola;
initv=[.5,2,2];
global alpha;
close all;
% alpha=0;
alpha=1e-4;
sv=['-o';'-*';'-v';'-s'];
hv=[1e-1,1e-2,1e-4,1e-8];
rand('state',1);

disp(sprintf('%s \t%s \t%s \t%s \t\t %s \t\t %s \t\t %s \t\t %s',...
	     'function','h','jdiff','jdiff2','x','f(x)','fp(x)','fpp(x)'));

for i=1:3,
  init=initv(i);
  f=sprintf('f%d',i);
  disp(sprintf(' '));
  for jdiff=0:1
    for jdiff2=jdiff:1      
      k=0;
      xh=zeros(20,4);

      %only do h loop if not pure analytic
      if jdiff+jdiff2 == 0
	hv2=0;
      else
	hv2=hv;
      end
      
      for h=hv2,
	[x, hist] = newtsol_opt(init, f, h, tola, tolr, jdiff, jdiff2);
	disp(sprintf('%s \t\t%g \t%d \t%d \t\t%12.5e \t%12.5e \t%12.5e \t%12.5e',...
		     f,h,jdiff,jdiff2,x,hist(end,2),hist(end,3),hist(end,4)));
	k=k+1;
	zh(1,k)=length(hist(:,1));
	xh(1:zh(1,k),k)=hist(:,1);
	yh(1:zh(1,k),k)=hist(:,3);
	wh(1:zh(1,k),k)=hist(:,2);
      end

      figure
      subplot(2,1,1)
      for l=1:k,
	semilogy(xh(1:zh(1,l),l),yh(1:zh(1,l),l),sv(l,:))
	hold on
      end
      if jdiff+jdiff2 ~= 0
	legend(sprintf('h=%g',hv2(1)),sprintf('h=%g',hv2(2))...
	       ,sprintf('h=%g',hv2(3)),sprintf('h=%g',hv2(4)));
      end
      
      hold off
      xlabel('Iteration number');
      ylabel('Gradient Norm');
      set(gca, 'LineWidth', 1.5);
      st=sprintf('fi=%d,jdiff=%d,jdiff2=%d,a=%g',i,jdiff,jdiff2,alpha);
      title(st);

      subplot(2,1,2)
      for l=1:k,
	semilogy(xh(1:zh(1,l),l),wh(1:zh(1,l),l),sv(l,:))
	hold on
      end
      if jdiff+jdiff2 ~= 0
	legend(sprintf('h=%g',hv2(1)),sprintf('h=%g',hv2(2))...
	       ,sprintf('h=%g',hv2(3)),sprintf('h=%g',hv2(4)));
      end
      
      hold off
      xlabel('Iteration number');
      ylabel('Function value');
      set(gca, 'LineWidth', 1.5);

      orient('tall');
      print('-deps2',sprintf('%s.eps',st));

      disp(sprintf(' '));
    end
  end
  
end

  
