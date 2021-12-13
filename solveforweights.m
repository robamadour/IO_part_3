
% This code takes the moments from the data, the standard deviation of the
% moments, and the value of the model moments at different parameter values
% and solves for the weights that set the weighted model moments closest to
% the data moments.  This uses "lsqlin" a least-squares linear solver with
% the constraints that the weights be between 0 and 1 and sum to 1.
 
% Written by Ashley, 12/18/2019

wd = "G:\Mi unidad\Columbia_doc\My Classes\IO I\Assignment\final_version\IO_partII";
cd(wd)

filename = 'Data\moment_data_class_assignment_use.csv';
alldata = csvread(filename,1);
% size(alldata)

filename2 = 'Data\moment_variance_class_assignment_use.csv';
delimiterIn = ',';
var = importdata(filename2,delimiterIn);


% for the assignment we're going to only use the variances, not the
% covariances. 
varsub = inv(var.*eye(size(var,1)));
%size(varsub)

sdmat = cholcov(varsub);
%size(sdmat)
 
datamoment = alldata(:,1);
%size(datamoment)

norm_data = sdmat*datamoment;
moments = alldata(:,2:end);
% size(moments)
norm_moment = sdmat*moments;
dimx = size(norm_moment,2);
 
A = ones(1,dimx);
b = 1;
Aeq = ones(1,dimx);
beq = 1;
lb = zeros(1,dimx);
up = ones(1,dimx); 
 
x0 = ones(dimx,1)./dimx;
 
tolx = 1*(10^(-10));
tolfun = 1*(10^-10);
evals = 1000000;
iters = 100000;
 
options = optimset('Display','iter','TolX',tolx,'TolFun',tolfun,'MaxFunEvals',evals,'MaxIter',iters);
[w,func,resid,exitflag,output,lambda] = lsqlin(norm_moment,norm_data,A,b,Aeq,beq,lb,up,x0,options );
 
running = cumsum(ones(size(w)));
 
a = running(w>0.0000001) ;
b = w(w>0.0000001);

out = [a b];

filenameout = 'Output\weightparams_class_assignment_P3Q5.csv';
csvwrite(filenameout,out);

%exit, clear
