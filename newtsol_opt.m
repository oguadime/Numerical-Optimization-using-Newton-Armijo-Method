function [x, hist] = newtsol_opt(x, f, h, tola, tolr, jdiff, jdiff2) 
%NEWTSOL_OBJ   Newton-Armijo code for minimizing scalar functions
%
% This code terminates on small relative-absolute errors
%
%         [X, HIST] = NEWTSOL(X, F, TOLA, TOLR, JDIFF, JDIFF2) 
%
% Inputs: X = initial iterate
%         F = function
%         TOLA = absolute error tolerance
%         TOLR = relative error tolerance
%         JDIFF = 0, analytic derivative 
%                    syntax: [f,fp]=f(x)
%         JDIFF = 1, forward difference derivative 
%         JDIFF2 same as for JDIFF but for second derivatives  
%
% Output: x=approximate result
%      HIST = array of iteration history, useful for tables and plots
%            The five columns are iteration number, function value, derivative, 
%            number of step size reductions done in the line
%            search, and the Newton iterate itself
%
%
%             If you leave out the hist array by calling newtsol as
%             Z = NEWTSOL(X, F) 
%                then storage for the history array is not allocated and
%                the iteration history is not stored.
%

if nargin < 3
  h=1e-8;
end
if nargin < 4
  tola = 1e-5;
end
if nargin < 5
  tolr = tola;
end
if nargin < 6
  jdiff = 1;
end
if nargin < 7
  jdiff2 = 1;
end

maxit=100;    % maximum iterations
maxitls=20;   % maximum iterations inside the line search
alpha=1.d-4; 
small=1.d-8;
large=1.d8;
%
% Initialize
%
itc=0;
xc=x;
[fc,fcp,fcpp]=fpeval(f,x,h,jdiff,jdiff2);
tol=tolr*abs(fcp)+tola;
%
% Store iteration history?
%

if nargout == 2
    histmp=zeros(maxit,6);
    histmp(itc+1,:)=[itc, fc, abs(fcp), fcpp, 0, x];
end
%
%    Main Newton loop
%
while(abs(fcp) > tol)
%
%    Iteration becoming unbounded?
%
    if(abs(x) > large*xc)
        disp(' iterate too large');
        if nargout == 2
            histmp(itc+1,:)=[itc, fc, abs(fcp), fcpp, 0, x];
            hist=histmp(1:itc+1,:);
        end
        return
    end
%
    lambda=1;
    iarm=0;           % line search iteration counter

%    Check for small relative derivatives.
    if(abs(fcpp) <= small*abs(fcp))
        disp(' small derivative error');
        if nargout == 2
            histmp(itc+1,:)=[itc, fc, abs(fcp), fcpp, 0, x];
            hist=histmp(1:itc+1,:);
	end
	return;
    end
    if(fcpp < 0)
        disp(' second derivative not positive ');
    end
     s=-fcp/fcpp;
    xt=x+lambda*s;
    [ft,ftp,ftpp]=fpeval(f,xt,h,jdiff,jdiff2);
%
%    Main line search loop. 
%
    while(abs(ftp) >= (1 - alpha*lambda)*abs(fcp) + 1.d-12)
        lambda=lambda/2;
        xt=x+lambda*s;
        [ft,ftp,ftpp]=fpeval(f,xt,h,jdiff,jdiff2);
        iarm=iarm+1;
%
%    Are you spending too much time in the line search?    
%
        if(iarm > maxitls)
            disp(' line search failure');
            if nargout == 2
                histmp(itc+1,:)=[itc, fc, abs(fcp), fcpp, 0, x];
                hist=histmp(1:itc+1,:);
            end
            return
        end
    end
%
%    Step accepted, continue the Newton iteration.
%
    x=xt;
    itc=itc+1;
%
%    Too many iterations?
%
    if(itc > maxit)
        disp(' maxit reached');
        if nargout == 2
            histmp(itc+1,:)=[itc, fc, abs(fcp), fcpp, 0, x];
            hist=histmp(1:itc+1,:);
        end
        return
    end
%
    [fc,fcp,fcpp]=fpeval(f,x,h,jdiff,jdiff2);
    if nargout == 2
       histmp(itc+1,:)=[itc, fc, abs(fcp), fcpp, iarm, x];
    end
end
%
%   Fix up the history array.
%
if nargout == 2
    hist=histmp(1:itc+1,:);
end

%   End of newtsol_opt.


%
%   Return the function value and its derivatives using finite
%   differences or an analytical result
%
function [fc,fp,fpp]=fpeval(f,x,h,jdiff,jdiff2)

  if nargin < 3
    h=1e-8;
  end
  if nargin < 4
    jdiff = 1;
  end
  if nargin < 5
    jdiff2 = 1;
  end
  
  if jdiff==1
    if jdiff2==1
      fc=feval(f,x);
      fpp=findiff2(f,x,fc,sqrt(h)); %numerical second derivative
    else
      [fc,fpp]=feval(f,x);
    end
    fp=findiff(f,x,fc,h); % use a numerical derivative
  else
    if jdiff2==1      
      [fc,fp]=feval(f,x);
      fpp=findiff2(f,x,fc,h); %numerical second derivative
    else
      [fc,fp,fpp]=feval(f,x);
    end
  end


%
% A simple forward difference in MATLAB.
%
function fp = findiff(f,x,fc,h)
  fright=feval(f,x+h);
  fp=(fright-fc)/h;

function fpp = findiff2(f,x,fc,h)
  fr=feval(f,x+h);
  fl=feval(f,x-h);
  fpp=(fr-2*fc+fl)/(h^2);

    
