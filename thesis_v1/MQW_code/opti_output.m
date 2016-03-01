function [stop,options,optchanged] = opti_output(options,optimvalues,flag)

persistent last_best;
persistent final_val;

optchanged = false;
stop = false;
switch flag
    case 'init' %output when started
        disp('Starting the algorithm');
        last_best = optimvalues.bestx;   
        final_val = optimvalues.bestfval;
        fprintf('Best Point:\n')
        fprintf('%.3f ',last_best);fprintf('\n');
        fprintf('Best Val:\n%.3df\n',final_val);
    case 'iter' % during the evolution
        best = optimvalues.bestfval;
        if(best<final_val) 
            last_best = optimvalues.bestx; 
            final_val =optimvalues.bestfval;
            fprintf('Best Point:\n')
            fprintf('%.3f ',last_best);fprintf('\n');
            fprintf('Best Val:\n%.3df\n',final_val);
        end
    case 'done'% output when finished
        disp('Performing final task');   
        fprintf('Best Point:\n')
        fprintf('%.3f ',last_best);fprintf('\n');
        fprintf('Best Val:\n%.3df\n',final_val);
    otherwise
        disp('Nothing to do');
end
end