function options = Settings_FIF2_v2(varargin)
% 
% SETTINGS_FIF2_V2 Constructs option structure for the algorithms
%
% EXAMPLES
%
% OPTIONS = SETTINGS_FIF2_V2 with no input arguments returns
%   setting structure with default values
%
% OPTIONS = SETTINGS_FIF2_V2('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
% OPTIONS = SETTINGS_FIF2_V2(OLDOPTIONS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTIONS.
%
%
% SETTINGS_FIF2_V2 PROPERTIES :
%
% GENERAL
% 
%   saveEnd          - (0) If value is 1 save outputs in "name file"_decomp_MIF_2D_vXX.mat (only at termination)
%   saveIntermediate - (0) If value is 1 save outputs in "name file"_inter_decomp_MIF_2D_vXX.mat
%                      every inner step
%
%   verbose          - (1) 0 the method is silent, 1 normal, >1 loud, 
%
%   plots            - (0) the algorithm does not produce plots,
%                      1 it produces them 
%   
%   saveplots        - (0) no saving
%                      
%
% SPECIFIC PARAMETERS
%
%  delta       (0.001) Stopping criterion
%  ExtPoints   (3)     Number of extrema allowed in the remainder
%  NIMFs       (1)     Number of IMFs we want to produce, not counting
%                         the remainder
%  alpha       ('ave') Parameter used for the mask length computation.
%                           Allowed values 'ave', [0,100] or 'Almost_min'.
%                           If set to 'ave' the mask length equals 
%                           round(2*Xi*(length of the signal)/(number of extrema)).
%                           If set to 0 the mask length is proportional to the
%                           0 percentile of the distances between two subsequent extrema.
%                           If set to 100 then it is proportional to the 
%                           100 percentile of the distances between two subsequent extrema. 
%                           Finally if set to 'Almost_min' it is set to 
%                           a value close to the minimum (30-th percentile).
%  Xi          (1.6)   Parameter we use to tune the mask length
%  extensionType  (based on wextend)
%                 'zpd'             Zero extension
%                 'sp0'	            Smooth extension of order 0
%                 'spd' or 'sp1'	Smooth extension of order 1
%                 'sym' or 'symh'	Symmetric padding (half point): boundary value symmetric replication
%                 'symw'	        Symmetric padding (whole point): boundary value symmetric replication
%                 'asym' or 'asymh'	Antisymmetric padding (half point): boundary value antisymmetric replication
%                 ('asymw')         Antisymmetric padding (whole point): boundary value antisymmetric replication
%                 'ppd'	            Periodized extension 
%                 'per'	            Periodized extension  - If the signal length 
%                                   is odd, wextend appends on the right a copy of the last value, and performs 
%                                   the extension using the 'ppd' mode. Otherwise, 'per' reduces to 'ppd'. 
%                                   This rule also applies to images.
%
%  MaxInner    (200)    Max number of inner steps
%  UseFFT      (true)   Boolean: if true, use the direct iteration by means of the FFT.
%  MonotoneMaskLength (true) Boolean: if true when the algorithm compute a new mask length that is smaller or equal 
%                            to the previous one then automatically it increases to 1.1 of the previous mask length. 
%                            If false it allows smaller mask lengths.
%
% ------------------------------------------------------
% EXAMPLE
%          
%   >> options = Settings_FIF2_v2('delta',0.08,'NIMFs',5,'plots',1) 
%   >> IMF = FIF2_v3(x,options)
%              
%  Executes algorithm FIF2_v3 with delta = 0.08, stop after we have at most 5 IMFs and a trend = 5, it produces plots.                            
% ------------------------------------------------------      
%
% See also FIF2_v3, WEXTEND
%
%
% Ref. Antonio Cicone, Haomin Zhou. "Multidimensional Iterative Filtering method
%      for the decomposition of high-dimensional non-stationary signals".
%      Cambridge Core in Numerical Mathematics: Theory, Methods and
%      Applications, Volume 10, Issue 2, Pages 278-298, 2017.
%      doi:10.4208/nmtma.2017.s05
%
%      Stefano Sfarra, Antonio Cicone, Bardia Yousefi, Stefano Perilli,
%      Leonardo Robol, Xavier P.V. Maldague.
%      "Maximizing the detection of thermal imprints in civil engineering
%      composites after a thermal stimulus - The contribution of an
%      innovative mathematical pre-processing tool: the 2D Fast Iterative
%      Filtering algorithm. Philosophy, comparisons, numerical, qualitative
%      and quantitative results". 2021. Submitted
% 

% (Ripped from sdpsettings.m by Johan Lufberg)


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help Settings_FIF2_v1
    return;
end


Names = {
    % General   
    'saveEnd'
    'saveIntermediate'
    'verbose'
    
    'maxTime'
    'plots'
    'saveplots'
     % 'saveinIt'
     % 'logfile'
     'algorithm'
     
    % FIF2
    'delta'
    'ExtPoints'
    'NIMFs'
    'Xi'
    'alpha'
    'extensionType'
    'MaxInner'
    'UseFFT'
    'MonotoneMaskLength'
};

obsoletenames ={ % Use when options have become obsolete
};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) && isstruct(varargin{1})
    options = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General 
    %options.saveinIt = 0;
    options.saveIntermediate = 0;
    options.saveEnd = 0;
    options.verbose = 1;
    %options.logfile = 1;
    options.maxTime = Inf;
    options.plots = 0.0;
    options.saveplots = 0; 
           
    % MIF
    options.delta = 0.001;
    options.ExtPoints=3;
    options.NIMFs=1;
    options.Xi=1.6;
    options.alpha='Almost_min';
    options.extensionType='asymw';
    options.MaxInner= 200;
    options.UseFFT=true;
    options.MonotoneMaskLength=true;
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;       % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability... No need yet
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;    % we expect a value next
    else
        eval(['options.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end

end

