classdef RoughnessStore
    methods (Static)
        function out = setgetIm(data)
            persistent im;
            if nargin 
                im = data;
            end
            out = im;
        end
    end
end