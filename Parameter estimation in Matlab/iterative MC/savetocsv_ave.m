
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

outf = [];
for Group = 3:4
    allPMs = alltres{Group,3};
    params = alltres{Group,2};
    gofset = alltres{Group,1};
    
    dims= size(allPMs);
    % for i = 1: dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                for l = 1:dims(4)
                    outf = [outf; Group j k l params(:,j,k,l)' gofset(j,k,l) allPMs(:,j,k,l)'];
                end
            end
        end
end
csvwrite('PMs_G34.csv',outf);

outf = [];
best = bestres;    
dims= size(best);
        for j = 1:dims(1)
            for k = 1:dims(2)
                for l = 1:dims(3)
                    outf = [outf; j k l squeeze(bestres(j,k,l,:))'];
                end
            end
        end

csvwrite('best_by_G34.csv',outf)