function [state, options,optchanged] = MyGAoutput(options,state,flag,p1,p2)
% customed output function to replace the one from GA toolbox 
% if the best individual is update, output the new individual and its fitness


% format short
persistent last_best;
persistent final_best;
persistent final_bestPop;
global para;
optchanged = false;

switch flag
    case 'init' %output when started
        disp('Starting the algorithm');
%         disp(state.Population(:,:));% print whole initial population
        % possible to insert customed initial individual here¡£
        [best,i] = min(state.Score);
        last_best = best;
        final_best = last_best;        
        final_bestPop = state.Population(i,:);
         
        output = round(decode(state.Population(i,:),para));
        str = sprintf('%-5.0f    %9.6e\n',...
               state.Generation,best);
        disp(str);
        disp(output');  
    case 'iter' % during the evolution
        [best,i] = min(state.Score);
        if last_best ~= best
            last_best = best;
            output = round(decode(state.Population(i,:),para));
            str = sprintf('%-5.0f    %9.6e\n',...
                   state.Generation,best);
            disp(str);
            disp(output');   
        end
        if(best<final_best) 
                final_best = best;
                final_bestPop = state.Population(i,:);
        end
    case 'done'% output when finished
        disp('Performing final task');   
        output = round(decode(final_bestPop,para));
		final_best
        disp([output']); 
    otherwise
        disp('Nothing to do');
end
% state = kicksimilar(state,3,10);% kick off some similar individuals.