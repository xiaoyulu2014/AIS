% adaptive nonparametric importance sampler
% Use normal estimate of average probability
% space is (0,1], binary split hierarchically, using ANIS within each split
% recursively
% Try to reduce KL loss

function RR = anis_1D(C)
% target density
%truep = @(x) (x>1).*(x-1) + .0;
    truep = @(x) exp(-10*(1-x));
    % threshold before splitting node
    threshold = 1;

    % binary hierarchy encoded as follows: root interval is at index 1, 
    % and for each interval at i, its two children are 2*i and 2*i+1, 
    % while its parent is floor(i/2).

    total = [truep(rand) 0 0]; % sum of estimates of mass in each interval
    dev = [total.^2 0 0];      % sum of square of estimates of mass in each interval
    nn = [1 0 0];              % number of samples in each interval
    mm = [total 0 0];          % current estimate of mass
    vv = [0 0 0];              % current variance of estimate
    explorefactpr = C;
    
    T = 10;           % number of proposals;
    cost = zeros(1,T);
    xx = zeros(1,T);     % samplesdev = [total.^2 0 0];      % sum of square of estimates of mass in each interval
     
    pp = zeros(1,T);     % p at samples
    maxlevel = 10;       % max depth of hierarchy (for computational reasons)
    RR = zeros(1,T);     % KL regret

        for tt=1:T
          % start at root
          ii = 1; % index of interval
          ll = 0; % level of hierarchy.  root is level 0
          ss = 0; % left-hand endpoint of interval we're sampling from

          qq = []; % sequence of proposal probabilities down levels
          % mm(ii) = total(ii)/nn(ii);
          % vv(ii) = (dev(ii)/nn(ii) - mm(ii)^2)/nn(ii);
          while(1)
            if (nn(2*ii)<threshold || nn(2*ii+1)<threshold) || ii>=2^maxlevel
              % don't recurse down hierarchy.  Uniform sample in interval
              rr = rand(1);
              xx(tt) = ss + rr/(2^ll);
              pp(tt) = truep(xx(tt));
              ii = 2*ii + (rr>.5);
              ll = ll + 1;
              qq(ll) = .5;
              break;
            else
              cc = [2*ii 2*ii+1]; % two children
              % construct proposal using hierarchy
              mm(cc) = total(cc)./nn(cc);
              vv(cc) = dev(cc)./nn(cc) - mm(cc).^2./nn(cc);
              % core step: construct exploration boosted proposal
              empiricalm = mm(cc);
              exploreboost = explorefactor*sqrt(log(nn(ii)+1)./(nn(cc)));
              rr = mm(cc) + exploreboost;
              % sample 
              ll = ll + 1;
              qq(ll) = rr(2) / sum(rr);
              jj = rand(1)<qq(ll);
              if jj==0, 
                qq(ll) = 1-qq(ll);
              end
              % recurse to next level
              ii = 2*ii+jj;
              ss = ss + jj/(2^ll);
              % if children nodes don't exist yet, make them.
              if 2*ii+1>length(nn)
                nn(2*ii+1) = 0;
                total(2*ii+1) = 0;
                dev(2*ii+1) = 0;
              end
            end
          end

          % at this point ii is index of node on tree, ll is level.
          % note that this is one level lower than level at which uniform proposal is
          % drawn, as we'd like to keep stats on the two children of the node as well.
          % xx(ii) is proposed point, pp(ii) is true density at xx(ii).
          % qq is sequence of proposal probabilities down levels

          rr = pp(tt);
          for kk = ll:-1:0
            % update suff stats
            total(ii) = total(ii) + rr;
            dev(ii) = dev(ii) + rr^2;
            nn(ii) = nn(ii) + 1;
            % up one level of tree
            if kk>0
              %rr = rr / qq(kk)/2;
              ii = floor(ii/2);
            end
          end

          % plot
        %  x = .01:.01:.99;
        %  p = truep(x);
        %  plot(x,p/.1,'r','linewidth',2);
        %  hold on
         tmp=[];
          yy = ones(1,length(nn));
          dd = zeros(1,length(nn));
          for ii = 1:length(nn)
            cc = [2*ii 2*ii+1];
            if 2*ii<=length(nn) && all(nn(cc)>threshold) && ii<2^maxlevel
              % construct proposal using hierarchy
              mm(cc) = total(cc)./nn(cc);
              rr = mm(cc) + explorefactor*sqrt(log(nn(ii)+1)./(nn(cc)));
              yy(cc) = yy(ii)*rr/sum(rr)*2;
              dd(ii) = 1;
            elseif ii==1||dd(floor(ii/2))==1
              yy(cc) = yy(ii);
              ll = floor(log2(ii));
              ss = (ii-2^ll)/2^ll;
         %     plot([ss ss+1/2^ll],[yy(ii) yy(ii)],'b','linewidth',2);
              YY = 0.1*(truep(ss+1/2^ll)-truep(ss));
              %hh = sum(xx(1:tt)>ss&xx(1:tt)<(ss+1/2^ll))/tt*2^ll;
              %plot([ss ss+1/2^ll],[hh hh],'g','linewidth',2);
              %regret
              %YY = zeros(1,2.^(maxlevel-ll));
              %for i= 1:(2.^(maxlevel-ll))
              %    YY(i) = 0.1*(exp(10*(ss+i/2^maxlevel-1))-exp(10*(ss+(i-1)/2^maxlevel-1)));
             % end
             % RR(tt) = RR(tt) + sum(YY.*log(YY./0.1./(repmat(yy(ii)/2.^(maxlevel-ll),1,2.^(maxlevel-ll))./2^ll)));
             RR(tt) = RR(tt)-log(yy(ii))*YY./0.1;
             tmp=[tmp,YY./0.1];
           %  tmp=[tmp,yy(ii)./2^ll];  this is  q_{at}
             %RR(tt) = RR(tt) + sum(YY./2.^ll.*log(YY./0.1./(yy(ii)./2^ll)));
            end
          end
          RR(tt) = RR(tt)+1.3030;
      %    hold off
      %    drawnow

        end
        
end



































