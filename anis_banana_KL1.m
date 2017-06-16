% adaptive nonparametric importance sampler
% Use normal estimate of average probability
% space is (0,1], binary split hierarchically, using ANIS within each split
% recursively
% Try to reduce KL loss
 
C=logspace(-2,2,30);

T = 1000;               % number of proposals
xx_mat=zeros(30,T);
xxy_mat=zeros(30,T);
RR_mat=zeros(30,T);


parfor j=1:30
% target density
%truep =  @(x1,x2) exp(-0.1*(0.01*x1^2+(x2+0.1*(x1^2-10))^2));
truep = @(x1,x2) exp(-0.5*(0.03*x2^2+(x1+0.03*(x2^2-100))^2));
% exploration factor
explorefactor = 5;
timer1 = 0;
timer2 = 0;
% binary hierarchy encoded as follows: root interval is at index 1, 
% and for each interval at i, its two children are 2*i and 2*i+1, 
% while its parent is floor(i/2).
M = 20;             % domain of X space
 
total = [truep(2*M*rand-1,2*M*rand-1) 0 0 ; 0 0 0 ; 0 0 0]; % sum of estimates of mass in each interval
nn = zeros(3,3);              % number of samples in each interval
mm = vertcat(horzcat(total,zeros(3,2)),zeros(2,5));          % current estimate of mass
dd = zeros(10000,10000);  % whether the node is a child

cost = zeros(1,T);
xx = [];     % samplesdev = [total.^2 0 0];      % sum of square of estimates of mass in each interval
xxy = [];     % samplesdev = [total.^2 0 0];      % sum of square of estimates of mass in each interval

pp = zeros(1,T);     % p at samples
maxlevel = 10;       % max depth of hierarchy (for computational reasons)
RR = zeros(1,T);     % KL regret
q = zeros(10000,10000); q(1,1)=1;

    for tt=1:T
      % start at root
      ii = 1; % index of interval
      iiy = 1; % index of interval
      ll = 0; % level of x-hierarchy.  root is level 0
      lly = 0; % level of y-hierarchy.  root is level 0
      ss = -M; % left-hand x-endpoint of rectangle we're sampling from
      ssy= -M;% lower-hand y-endpoint of rectangle we're sampling from

      qq = []; % sequence of proposal probabilities down levels
      qqy = []; % sequence of proposal probabilities down levels
      % mm(ii) = total(ii)/nn(ii);
      % vv(ii) = (dev(ii)/nn(ii) - mm(ii)^2)/nn(ii);
      while(1)   
          if dd(ii,iiy) == 0 
             tic;
             rr = rand(1); rry = rand(1); 
             xx(tt) = ss + 2*M*rr/(2^ll);
             xxy(tt) = ssy + 2*M*rry/(2^lly);
             timer1 = (toc + (tt-1)*timer1) / tt;
             pp(tt) = truep(xx(tt),xxy(tt));
             
            if ll==lly
                
                indexx1 = find(xx>ss & xx < ss+M/(2^ll));
                indexx2 = find(xx > ss+M/(2^ll) & xx < ss+2*M/(2^ll));
                indexy = find(xxy> ssy & xxy < ssy+2*M/(2^lly));
                indexx = intersect(union(indexx1,indexx2),indexy);

               index1=intersect(indexx1,indexy); 
               index2=intersect(indexx2,indexy);
               rr=[sum(pp(index1)) sum(pp(index2))]./[length(index1) length(index2)]  ...
                        + explorefactor*sqrt(log(nn(ii,iiy)+1))./[length(index1) length(index2)];
               q_tmp = rr/sum(rr);   
               tic; q_tmp > rand(1);
               timer2=(toc + (tt-1)*timer2) / tt;
              % weights1 = pp(indexx)/(q(ii,iiy)*4*M^2/2^(ll+lly));
               if isempty(index1)==0 & isempty(index2)==0 & nn(ii,iiy)>100
                  % weights21 =  pp(index1)./(q_tmp(1)*8*M^2/2^(ll+lly));
                  % weights22 =  pp(index2)./(q_tmp(2)*8*M^2/2^(ll+lly));
                   A0 = (2*M)^2/2^(ll+lly);
                   %KL1 = sum(pp(indexx).*log(A0*pp(indexx)))/sum(pp(indexx));
                   %KL2 = sum(pp(index1).*log(A0*pp(index1)/(2*q_tmp(1))))/sum(pp(indexx))...
                     %+sum(pp(index2).*log(A0*pp(index2)/(2*q_tmp(2))))/sum(pp(indexx)); 
                 %KLdiff = KL2-KL1;
                 var_tmp = var([pp(index1)*log(2*q_tmp(1))./mean(pp(indexx)) pp(index2)*log(2*q_tmp(2))./mean(pp(indexx))]);
                 mean_tmp = mean([pp(index1)*log(2*q_tmp(1))./mean(pp(indexx)) pp(index2)*log(2*q_tmp(2))./mean(pp(indexx))]);
                 stats = (mean_tmp + log(timer1/(timer1+timer2)))/sqrt(var_tmp);
                 %stats=length(find([pp(index1)*log(2*q_tmp(1)) pp(index2)*log(2*q_tmp(2))] > -mean(pp(indexx))*log(timer1/(timer1+timer2))));
                   if  stats > 3.19 %stats/length(indexx) > 0.005 %reject null %KLdiff < log(timer1/(timer1+timer2))
                       dd(ii,iiy) = 1; %split, (ii,iiy) becomes parent node;
                       total([2*ii 2*ii+1],iiy) =  [sum(pp(index1)) sum(pp(index2))]*2*M^2/2^(ll+lly);
                       nn([2*ii 2*ii+1], iiy) = [length(index1) length(index2)];
                   end
               end       
                   
                % update statistics
               kk=ll;kky=lly;
               while kk>=0
                  while kky>=0
                    total(ii,iiy) = total(ii,iiy) + pp(tt)*(2*M/2^kk)*(2*M/2^kky);
                    nn(ii,iiy) = nn(ii,iiy) + 1;
                    if (kky==kk & kky>0)
                     iiy = floor(iiy/2);
                     kky=kky-1;
                    elseif (kky<kk & kk>0)
                     ii = floor(ii/2);
                     kk=kk-1;
                    end
                    if kk==0 & kky==0
                        nn(1,1)=nn(1,1)+1;
                    break
                    end
                  end
                  if kk==0 & kky==0
                        nn(1,1)=nn(1,1)+1;
                    break
                  end
                end
               break
          elseif ll>lly
                indexx = find(xx>ss & xx < ss+2*M/(2^ll));
                indexy1 = find(xxy >ssy & xxy< ssy+M/(2^lly));
                indexy2 = find(xxy> ssy+M/(2^lly) & xxy < ssy+2*M/(2^lly));
                indexx = intersect(union(indexy1,indexy2),indexx);

               index1=intersect(indexy1,indexx);
               index2=intersect(indexy2,indexx);
               rr=[sum(pp(index1)) sum(pp(index2))]./[length(index1) length(index2)]  ...
                        + explorefactor*sqrt(log(nn(ii,iiy)+1))./[length(index1) length(index2)];
               q_tmp = rr/sum(rr);   
               tic; rry>0.5;
               timer2=(toc + (tt-1)*timer2) / tt;
               %weights1 = pp(indexx)/(q(ii,iiy)*4*M^2/2^(ll+lly));
               if isempty(index1)==0 && isempty(index2)==0  & nn(ii,iiy)>100
                  % weights21 =  pp(index1)./(q_tmp(1)*8*M^2/2^(ll+lly));
                  % weights22 =  pp(index2)./(q_tmp(2)*8*M^2/2^(ll+lly));
                 A0 = (2*M)^2/2^(ll+lly);
                   %KL1 = sum(pp(indexx).*log(A0*pp(indexx)))/sum(pp(indexx));
                   %KL2 = sum(pp(index1).*log(A0*pp(index1)/(2*q_tmp(1))))/sum(pp(indexx))...
                     %+sum(pp(index2).*log(A0*pp(index2)/(2*q_tmp(2))))/sum(pp(indexx)); 
                 %KLdiff = KL2-KL1;
                 var_tmp = var([pp(index1)*log(2*q_tmp(1))./mean(pp(indexx)) pp(index2)*log(2*q_tmp(2))./mean(pp(indexx))]);
                 mean_tmp = mean([pp(index1)*log(2*q_tmp(1))./mean(pp(indexx)) pp(index2)*log(2*q_tmp(2))./mean(pp(indexx))]);
                 stats = (mean_tmp + log(timer1/(timer1+timer2)))/var_tmp;
              %stats=length(find([pp(index1)*log(2*q_tmp(1)) pp(index2)*log(2*q_tmp(2))] > -mean(pp(indexx))*log(timer1/(timer1+timer2))));
                   if  stats > 3.19 %stats/length(indexx) > 0.005 %reject null %KLdiff < log(timer1/(timer1+timer2))
                       dd(ii,iiy) = 1; %split, (ii,iiy) becomes parent node;
                       total(ii,[2*iiy 2*iiy+1]) =  [sum(pp(index1)) sum(pp(index2))]*2*M^2/2^(ll+lly);
                       nn(ii,[2*iiy 2*iiy+1]) = [length(index1) length(index2)];
                   end
               end        
                % update statistics
               kk=ll;kky=lly;
               while kk>=0
                  while kky>=0
                    total(ii,iiy) = total(ii,iiy) + pp(tt)*(2*M/2^kk)*(2*M/2^kky);
                    nn(ii,iiy) = nn(ii,iiy) + 1;
                    if (kky==kk & kky>0)
                     iiy = floor(iiy/2);
                     kky=kky-1;
                    elseif (kky<kk & kk>0)
                     ii = floor(ii/2);
                     kk=kk-1;
                    end
                    if kk==0 & kky==0
                        nn(1,1) = nn(1,1)+1;
                    break
                    end
                  end
                  if kk==0 & kky==0
                      nn(1,1) = nn(1,1)+1;
                    break
                    end
                end
               break
            end
        else
          if ll==lly+1
          % recurse to next y-level
              ccy = [2*iiy 2*iiy+1]; % two children
              % construct proposal using hierarchy
              mm(ii,ccy) = total(ii,ccy)./(nn(ii,ccy));
              % core step: construct exploration boosted proposal
              exploreboost = explorefactor*sqrt(log(nn(ii,iiy)+1))./nn(ii,ccy);
              rry = mm(ii,ccy) + exploreboost;
              % sample 
              lly = lly + 1;
              qqy(lly) = rry(2) / sum(rry);
              q(ii,ccy) = rry/sum(rry)*q(ii,iiy);
              jjy = rand(1)<qqy(lly);
              if jjy==0, 
                qq(lly) = 1-qqy(ll);
              end
              % recurse to next level
              iiy = 2*iiy+jjy;
              ssy = ssy + jjy/(2^ll)*2*M;
              % if children nodes don't exist yet, make them.
              if 2*iiy+1>size(nn,2)
                nn(:,2*iiy+1) = 0;
                total(:,2*iiy+1) = 0;
              end
          elseif  ll==lly
              cc = [2*ii 2*ii+1]; % two children
              % construct proposal using hierarchy
              mm(cc,iiy) = total(cc,iiy)./(nn(cc,iiy));
              % core step: construct exploration boosted proposal
              exploreboost = explorefactor*sqrt(log(nn(ii,iiy)+1))./nn(cc,iiy);
              rr = mm(cc,iiy) + exploreboost;
              % sample 
              ll = ll + 1;
              qq(ll) = rr(2) / sum(rr);
              q(cc,iiy) = rr/sum(rr)*q(ii,iiy);
              jj = rand(1)<qq(ll);
              if jj==0, 
                qq(ll) = 1-qq(ll);
              end
              % recurse to next level
              ii = 2*ii+jj;
              ss = ss + jj/(2^ll)*2*M;
              if 2*ii+1>size(nn,1)
                nn(2*ii+1,:) = 0;
                total(2*ii+1,:) = 0;
              end
          end
        end
      end
    


      % at this point ii is index of node on tree, ll is level.
      % note that this is one level lower than level at which uniform proposal is
      % drawn, as we'd like to keep stats on the two children of the node as well.
      % xx(ii) is proposed point, pp(ii) is true density at xx(ii).
      % qq is sequence of proposal probabilities down levels
     % if 0
        %plot
     %  ;     
     if rem(tt,10000)==0
         try
            subplot(3,3,tt/1000)
            %figure()
            yy=ones(size(nn));
            for ii=1:size(nn,1)
                for jj=1:size(nn,2)
                     cc = [2*ii 2*ii+1];  ccy = [2*jj 2*jj+1];
                     ll = floor(log2(ii));lly = floor(log2(jj));
                     if  ll==lly &&  2*ii < size(nn,1)  && dd(ii,jj)==1 && ii < 2^maxlevel
                         rr = total(cc,jj)./(nn(cc,jj)) + explorefactor*sqrt(log(nn(ii,jj)+1))./nn(cc,jj);
                         yy(cc,jj)=rr./sum(rr).*yy(ii,jj);
                     elseif ll==lly+1 && 2*jj <size(nn,2) && dd(ii,jj)==1  && jj < 2^maxlevel
                         rr = total(ii,ccy)./nn(ii,ccy) + explorefactor*sqrt(log(nn(ii,jj)+1))./nn(ii,ccy);
                         yy(ii,ccy)=rr./sum(rr).*yy(ii,jj);
                     elseif  nn(ii,jj)>0 && dd(ii,jj)==0 && ((ll==lly && dd(ii,floor(jj/2))==1) || (ll==lly+1 && dd(floor(ii/2),jj)==1))
                         ss = (ii-2^ll)/2^ll*(2*M)-M; ssy = (jj-2^lly)/2^lly*(2*M)-M;
                       %  surf([ss,ss+2*M/(2^ll)],[ssy,ssy+2*M/(2^lly)],yy(ii,jj)/(2*M/(2^ll)*2*M/(2^lly))*ones(2,2),'edgecolor','g');hold on; 
                         surf([ss,ss+2*M/(2^ll)],[ssy,ssy+2*M/(2^lly)],yy(ii,jj)/(2*M/(2^ll)*2*M/(2^lly))*ones(2,2),'edgecolor',[0 0.7 0]);view(2); hold on;  
                         x_mid = ss + M/(2^ll); y_mid = ssy + M/(2^lly);                 
                         RR(tt) = RR(tt) + (truep(x_mid,y_mid)*log(truep(x_mid,y_mid))...
                             - log(yy(ii,jj)/(2*M/(2^ll)*2*M/(2^lly)))) * 2*M/(2^ll)*2*M/(2^lly);

                     end
                end
            end
            colormap summer
            hold off
            drawnow 
            colorbar()
         end
     end
    end

     
    xx_mat(j,:)=xx;
    xxy_mat(j,:)=xxy;
    RR_mat(j,:) = RR;

    end

    
figure();hist3(horzcat(xx',xxy'),[50,50])




save('banana/result.mat')
    
       
         
%sqrt(log(t))/n  best C(12)







y = -20:1:20; x= -20:1:20;
for i=1:length(x)
  for j=1:length(y)
      p(i,j) = truep(x(i),y(j));
  end
end
figure();surf(x,y,p');colorbar;%view(2);
xlabel('x');ylabel('y');


for i=1:30
    figure();hist3(horzcat(xx_mat(i,:)',xxy_mat(i,:)'),[100,100])
    xlabel('x');ylabel('y')
end
    
for i=1:30
    plot(RR_mat(i,:));hold on;
end

  [X,Y]=meshgrid([ss,ss+2*M/(2^ll)],[ssy,ssy+2*M/(2^lly)]);
  figure();surf(X,Y,yy(ii,jj)/(2*M/(2^ll)*2*M/(2^lly))*ones(2,2),'edgecolor','none'); hold on; view(2); 


  
  
  
for j=1:10
    plot((regret(j,:)));hold on
end
legend('1','2','3,','4','5','6','7','8','9','10','location','southeast')

plot(mean(regret,1))



R1 = mean(regret,1)
plot(cumsum(R1(400:end)),'r');hold on ;
plot(cumsum(R2(400:end)));hold on ;
plot(cumsum(R3(400:end)));hold on ;
plot(cumsum(R4(400:end)));hold on ;
plot(cumsum(R5(400:end)));hold on ;
plot(cumsum(R6(400:end)));hold on ;

legend('sqrt(log(t)/n)','sqrt(log(t))/n','sqrt(log(t))/n^{0.25}','log(t)/sqrt(n)','t^{0.25}/n^{0.5}',...
'log(t)/n','Location','southeast')
title('KL cumulative regret')
xlabel('number of iteration')
ylabel('regret')




boost=['sqrt(log(t)/n)','sqrt(log(t))/n','sqrt(log(t))/n^0.25','log(t)/sqrt(n)','t^0.25/n^0.5']




save('regret_hierarchy.mat','R1','R2','R3','R4','R5','R6')























