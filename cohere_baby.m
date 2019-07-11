%input is data in the time domain

function [outmat, F] = cohere_baby(inmat);

for chan1 = 1:size(inmat,1)
    disp('running channel: ')
    disp(chan1)
   for chan2 = 1:size(inmat,1)
       Cxy = []; 
     for segment = 1:size(inmat,3)
       
        X = squeeze(inmat(chan1, :, segment));
        Y = squeeze(inmat(chan2, :, segment));
        
        [Cxy,F] = mscohere(X,Y, 100, [], 500, 500); % this the whole frequency spectrum as coherence vals
        
        if segment == 1
            sumCxy = Cxy; 
        else
            sumCxy = sumCxy+Cxy; 
        end
        
     end
     outmat(chan1, chan2, :) = sumCxy./segment;
   end

end
  outmat = outmat(:,:,1:21);
  F = F(1:21);
end
