#!/home/wwalten/bin/octave -q
pre="lastrun";
if (exist ("residua")) 
	pre=residua;
end
if (nargin==1)
	pre()=nth(argv,1);
end
disp(["Loading " pre ".*"]);
fits=["riemann";"global";"kalman"];
clear detcov;
ind=0; detcov=[];
for fit=1:3
	for num=0:3
		filename=[pre "." deblank(fits(fit,:)) int2str(num)];
		[op stat]=system(["ls " filename " >& /dev/null" ]);
		if (~stat)
			ind=ind+1;
			eval (["clear " pre]);
			load (filename);
%			disp(["loading " filename]);
			eval([ "res_" deblank(fits(fit,:)) int2str(num) " =" pre ]);
			lastrun=eval(pre);
			tmp=sprintf("det(cov(%.6s%d))= %e",deblank(fits(fit,:)),num,det(cov(lastrun)))
			tmp=sprintf("std(%.6s%d)= %e %e %e %e %e",deblank(fits(fit,:)),num,std(lastrun))
			clear lastrun;
			detcov=[detcov; tmp ];
		end
	end
end

eval (["clear " pre]);
clear num pre fits fit filename op stat ind tmp;
disp('done');
