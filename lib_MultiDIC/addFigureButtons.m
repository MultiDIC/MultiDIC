function []=addFigureButtons(varargin)

nargin=numel(varargin);

switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
    otherwise
        error('wrong number of input arguments');
end

addColorbarLimitsButton(hf);
addColormapButton(hf);
addEdgeColorButton(hf);
addFaceAlphaButton(hf);
addLightButton(hf);
addAmbientStrengthButton(hf);
addDiffuseStrengthButton(hf);
addSpecularStrengthButton(hf);
addFaceLightingButton(hf);

end