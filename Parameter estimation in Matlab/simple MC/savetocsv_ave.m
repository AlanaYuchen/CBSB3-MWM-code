
% outf = [];
% for mouse = 1:12
%     gofset= alltres{mouse,1};
%     params = alltres{mouse,2};
%     
%     dims= size(gofset);
%     for i = 1: dims(1)
%         for j = 1:dims(2)
%             for k = 1:dims(3)
%                outf = [outf; mouse i j k params(:,i,j,k)' gofset(i,j,k)];
%             end
%         end
%     end
% end
% csvwrite('S1_PC1.csv',outf);


% output PMs
%
%outf = [];
%for Group = 1:4
%    allPMs = alltres{Group,3};
%    params = alltres{Group,2};
%   gofset = alltres{Group,1};
    
%    dims= size(allPMs);
    % for i = 1: dims(1)
%        for j = 1:dims(2) %npset
%            for k = 1:dims(3) %day
%                for l = 1:dims(4) %trial
%                    outf = [outf; Group j k l params(:,j,k,l)' gofset(j,k,l) allPMs(:,j,k,l)'];
%                end
%           end
%       end
%end
%csvwrite('PMs_simpleMC.csv',outf);

outf = [];
best = bestres;    
dims= size(best);
        for j = 1:4 %group
            for k = 1:dims(2) %day
                for l = 1:dims(3) %trial
                    outf = [outf; j k l squeeze(bestres(j,k,l,:))'];
                end
            end
        end

csvwrite('best_S2T12.csv',outf)