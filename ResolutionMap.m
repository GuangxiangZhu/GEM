function [ ind ] = ResolutionMap( res )

% Map resolutions into bins (10k,50k,100k,150k,500k,1M). 
resolution=[10E3,50E3,100E3,250E3,500E3,1E6];
ind=size(resolution,2);
for i=1:size(resolution,2)-1
    if res<=(resolution(i)+resolution(i+1))/2
        ind=i;
        break;
    end 
end

end

