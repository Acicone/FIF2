function [IMF,logM] = FIF2_v3(f,options,M)

%
% Generate the decomposition of a nD signal f :
%
%  f = IMF(:,:,1) + IMF(:,:,2) + ... + IMF(:, :,size(IMF, 3))
%
% where the last component is the trend and other components are IMFs
% The mask length is computed based on alpha and the mask length values
% from the horizontal middle section of the signal
%
%                    Input
%
%   f          nD signal to be decomposed
%
%   options    Structure, generated using function Settings_FIF2_v2, containing
%              all the parameters needed in the various algorithms
%
%
%   See also SETTINGS_FIF2_V2, MAXMINS_V3_8, GETMASK_2D_V3.
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


%% deal with the input

if nargin == 0,  help FIF2_v3; IMF=[];logM=[];return; end
if nargin == 1, options = Settings_FIF2_v2; end
if nargin == 2, M=[]; end
%FigCol = 'ckmygr'; % Plot Colors
tol=10^-12;
% we extend the signal outside the boundaries to reduce the errors at the
% boundaries

N_o = size(f);
if strcmp(options.extensionType,'per')
    Ni=0;
else
    Ni=round(max(N_o)/2);
    
    f = wextend('2D',options.extensionType,f,Ni);
end
N = size(f);


IMF =[];

if options.saveplots>0
    nameFile=input('Please enter the name of the file as a string using '' and '' <<  ');
else
    nameFile=['Test_' datestr(clock,'YYYY_mm_DD_HH_MM_SS')];
end

%% Main code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Iterative Filtering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('prefixed_double_filter','MM');

nIMFs=0;
logM=[];
while nIMFs < options.NIMFs
    nIMFs=nIMFs+1;
    if options.verbose>0
        fprintf('\n IMF # %1.0d\n',nIMFs)
    end
    SD=1;
    SDlog=[];
    havelog=[];
    h=f;
    
    %%%%%%%%%%%%%%% Identify extreme points of the vectorized signal %%%%%%%%%%%%%%
    maxmins_f=Maxmins_v3_8(f(:),tol);
    %         figure
    %         plot(f(round(end/2),:))
    diffMaxmins_f=diff(maxmins_f);
    k = length(maxmins_f);
    
    
    
    
    if length(M)<nIMFs
        if k<=options.ExtPoints
            disp(['Number of Extrema smaller than ExtPoints = ' num2str(options.ExtPoints)])
            IMF(:,:,end+1) = f;
            IMF=IMF(Ni+1:end-Ni,Ni+1:end-Ni,:);
            return
        end
        if isa(options.alpha,'char')
            if strcmp(options.alpha,'ave') % Using an average mask length
                m = round(2*(length(f(:))/k)*options.Xi);
            elseif strcmp(options.alpha,'Almost_min') % Using an almost min mask length
                if 2*round(options.Xi*prctile(diffMaxmins_f,30))<2*round(length(f(:))/k*options.Xi)
                    m = 2*round(options.Xi*prctile(diffMaxmins_f,30));
                else
                    m = 2*round(length(f(:))/k*options.Xi);
                end
            elseif strcmp(options.alpha,'Median') % Using a median mask length                
                    m = 2*round(options.Xi*median(diffMaxmins_f));      
            else
                disp(' Value of alpha not recognized')
                IMF=[];
                logM=[];
                return
            end
        else % using a prefixed value alpha
            m = 2*round(options.Xi*prctile(diffMaxmins_f,options.alpha));
        end
        if nIMFs>1
            if m<=logM(nIMFs-1)
                fprintf('Warning mask length is decreasing at step %1d. ',nIMFs)
                if options.MonotoneMaskLength==true
                    m=ceil(logM(nIMFs-1)*1.1);
                    fprintf('The old mask length is %1d whereas the new one is forced to be %1d.\n',logM(nIMFs-1),ceil(logM(nIMFs-1)*1.1))
                else
                    fprintf('The old mask length is %1d whereas the new one is %1d.\n',logM(nIMFs-1),m)
                end
            end
        end
    else
        m=M(nIMFs);
    end
    
    %         fprintf('Value of the mask length = %1.0d\n',m)
    %         m=input('Input a new value \n<< ');
    
    
    logM(nIMFs)=m;
    
    
    
    inStepN=0;
    if options.verbose>0
        fprintf('\n  step #            SD             Mask length \n\n')
    end
    
    A = get_mask_2D_v3(MM,m);
    if options.plots>0
        ImgIMF=figure;
        ImgR=figure;
    end
    
    if options.UseFFT
        % Precompute the FFT of the signal and of the filter, so we can
        % apply it with almost no computations
        [fftm, fftn] = size(h);
        [fftp, fftq] = size(A);
        
        % Note that if the filter has a support larger than our image
        % we might need to extend the image a little bit. The repmat
        % command takes care of that.
        if fftp > fftm || fftq > fftn
            fftExtL=ceil(max(fftp/fftm, fftn/fftq));
            fftH = repmat(h,fftExtL);
            fftm_o=fftm;
            fftn_o=fftn;
            [fftm, fftn] = size(fftH);
        else
            fftm_o=fftm;
            fftn_o=fftn;
            fftH = h;
        end
        
        % Pad A the right way -- this is required to make sure the
        % center of gravity of the filter is in position (1,1) before
        % taking the FFT
        l1 = floor(fftp/2); m1 = floor(fftq/2);
        l2 = fftp - l1; m2 = fftq - m1;
        fftA = zeros(fftm, fftn);
        fftA(1:l2,1:m2) = A(l1+1:end,m1+1:end);
        fftA(end-l1+1:end,end-m1+1:end) = A(1:l1,1:m1);
        fftA(end-l1+1:end,1:m2) = A(1:l1,m1+1:end);
        fftA(1:l2,end-m1+1:end) = A(l1+1:end,1:m1);
        
        % Precomputing FFTs for later use
        fftH = fft2(fftH);
        fftA = fft2(fftA);
        
        % r1 and r2 are updated throughout the iterations so that r1 /
        % r2 is equal to the relative change in the solution at step j.
        % To accomplish these, some values are accumulated in
        % filter_factors (which indeed describes the action of the
        % accumulated filter on all the frequencies), making use of the
        % variable incr_filter.
        r2 = abs(fftH).^2;
        r1 = r2 .* abs(fftA).^2;
        incr_filter = (1 - abs(fftA)).^2;
        filter_factors = ones(size(fftA));
    end
    
    
    while SD>options.delta && inStepN < options.MaxInner
        inStepN=inStepN+1;
        
        
        % FIXME: Similar to the 1D case, this might be made more
        % efficient skipping the iterations completely, and just
        % computing the necessary index j.
        if options.UseFFT
            % Construct the residual for checking the stopping
            % criterion
            h_avenrm = sqrt(sum(sum( r2 .* filter_factors )));
            SD = sum(sum( r1 .* filter_factors )) ./ ...
                h_avenrm^2;
            
            % We never explicitly compute h_ave, and we construct the
            % final solution h only at convergence
            if SD < options.delta || inStepN >= options.MaxInner
                h = real(ifft2( fftH .* (1 - fftA).^(inStepN) ));
            else
                % Update the filter factors for the next round
                filter_factors = incr_filter .* filter_factors;
            end
        else
            h_ave = conv2(h, A, 'same');
            
            h_avenrm = norm(h_ave, 'fro');
            SD=h_avenrm^2 / norm(h, 'fro')^2;
        end
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        
        SDlog=[SDlog SD];
        havelog=[havelog h_avenrm];
        if options.verbose>0
            fprintf('    %2.0d      %1.14f          %2.0d\n',inStepN,SD,m)
        end
        
        %%%%%%%%%%%%%%%%%% generating f_n %%%%%%%%%%%%%%%%%%
        
        % If using the FFT, then we do not need to subtract the
        % average, which is not computed explicitly. The function h
        % will be automatically constructed at convergence by a direct
        % approach applying the filter with FFT (1 - fftA).^j.
        if ~options.UseFFT
            h=h-h_ave;
        end
        
        if options.plots>0
            figure(ImgIMF)
            surfH=surf(h);
            set(surfH, 'edgecolor','none')
            title('IMF')
            pause(0.5)
            %%
            figure(ImgR)
            surfF=surf(f-h);
            set(surfF, 'edgecolor','none')
            title('Remainder')
            pause(0.5)
        end
        if options.saveIntermediate == 1
            save([ nameFile '_inter_FIF2_v3.mat'])
        end
        
        if inStepN >= options.MaxInner && options.verbose>0
            disp('Max # of inner steps reached')
        end
    end
    
    
    
    %close all
    
    IMF(:,:,nIMFs) = h(1:fftm_o,1:fftn_o);
    f=f-h(1:fftm_o,1:fftn_o);
    
    %         figure
    %         plot(f(round(end/2),:))
    %         pause
end

IMF(:,:,end+1) = f;


% we remove the initial extension
IMF=IMF(Ni+1:end-Ni,Ni+1:end-Ni,:);

if options.saveEnd == 1
    save([ nameFile '_FIF2_v3.mat'])
end

end


%% Auxiliar functions


function A=get_mask_2D_v3(w,k)
%
%% get the mask with length 2*k+1 x 2*k+1
% k must be integer
% w is the area under the curve for each bar
% A  the mask with length 2*k+1 x 2*k+1

L=length(w);
m=(L-1)/2;  %2*m+1 =L punti nel filtro
w=[w zeros(1,(L-1)/2)];
A=zeros(k+1,k+1);
if k<=m % The prefixed filter contains enough points
    if mod(k,1)==0     % if the mask_length is an integer
        for i=0:k
            for j=0:k
                d1=sqrt(i^2+j^2);
                d=k;
                if d1 > d
                    A(j+1,i+1)=0;
                else
                    s=(m-1)+L/2*d1/d;
                    t=s+2;
                    %             N=(L-1)/(2*d);
                    %             s=(m+1)+(d1-1/2)*N;
                    %             t=(m+1)+(d1+1/2)*N;
                    s2=ceil(s)-s;
                    t1=t-floor(t);
                    A(j+1,i+1)=sum(w(ceil(s):floor(t)))+s2*w(ceil(s))+t1*w(floor(t));
                end
            end
        end
        A=[rot90(A(2:end,2:end),2) rot90(A(2:end,:),3)';rot90(A(:,2:end),1)' A];
        A=A/sum(sum(A,1),2);
    else   % if the mask length is not an integer
        disp('Need to write the code!')
        A=[];
        return
    end
else % We need a filter with more points than MM, we use interpolation
    disp('Need to write the code!')
    A=[];
    return
end
end





function varargout = Maxmins_v3_8(f,tol)

% v3_7 Based on version 3_6
% Minor revision: we extend f as df in order to make possible to compute
% the relative derivative like df(i)/abs(f(i)) for every i=1:N
%
% v3_6 Based on version 3_5
% Minor revision: we consider relative differences to avoid having problems with small or big values in a signal
%
% v3_5 Based on version 3_4
% Minor revision: we consider only periodical extension
%
% v3_4 Based on version 3
% Minor revisions: 1) added for constant extention the checking for Mins and
%                     Maxs emptiness
%                  2) completed the code for the periodical case
%
% v3 is Based on Version 2.
% Modified the way zero-derivative regions are handled.
%
% Identify the maxima and minima of a signal f

N = length(f);
Maxs = zeros(1,N);
Mins = zeros(1,N);
df = diff(f);

if size(f,1)>size(f,2)
    f=f.';
end

h = 1;

while h<N && abs(df(h)/f(h)) <= tol
    h=h+1;
end
if h==N
    if nargout<=1
        varargout{1}=[];
    elseif nargout==2
        varargout{1}=[];
        varargout{2}=[];
    end
    return
end

cmaxs=0;
cmins=0;

c = 0;

N_old=N;

df=diff([f f(2:h+1)]);
f=[f f(2:h)];
N=N+h;


last_df=[];
for i=h:N-2
    if   df(i)*df(i+1)/abs(f(i))^2 <= tol && df(i)*df(i+1)/abs(f(i))^2 >= -tol
        if df(i)/abs(f(i)) < -tol
            last_df=-1;
            posc = i;
        elseif df(i)/abs(f(i)) > tol
            last_df=+1;
            posc = i;
        elseif df(i)==0
            last_df=0;
            posc = i;
        end
        c = c + 1;
        if df(i+1)/abs(f(i)) < -tol
            if last_df==+1 || last_df==0
                cmaxs=cmaxs+1;                
                Maxs(cmaxs)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        if df(i+1)/abs(f(i)) > tol
            if last_df==-1 || last_df==0
                cmins=cmins+1;
                Mins(cmins)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        
    end
    if   df(i)*df(i+1)/abs(f(i))^2 < -tol
        if df(i)/abs(f(i)) < -tol && df(i+1)/abs(f(i)) > tol
            cmins=cmins+1;
            Mins(cmins)=mod(i+1,N_old);
            if Mins(cmins)==0
                Mins(cmins)=1;
            end
            last_df=-1;
        elseif df(i)/abs(f(i)) > tol && df(i+1)/abs(f(i)) < -tol
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=mod(i+1,N_old);
            if Maxs(cmaxs)==0
                Maxs(cmaxs)=1;
            end
            last_df=+1;
        end
    end
end
if c > 0
    %         % we deal with the boundary
    %         df_0=f(N)-f(1);
    %         if df_0==0
    %             if Initial_df < 0
    %                 if last_df==+1
    %                     cmaxs=cmaxs+1;
    %                     Maxs(cmaxs)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             elseif Initial_df > 0
    %                 if last_df==-1
    %                     cmins=cmins+1;
    %                     Mins(cmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             end
    %         else
    %             disp('Code missing!')
    %         end
    if cmins>0 && Mins(cmins)==0
        Mins(cmins)=N;
    end
    if cmaxs>0 && Maxs(cmaxs)==0
        Maxs(cmaxs)=N;
    end
end

Maxs=Maxs(1:cmaxs);
Mins=Mins(1:cmins);
maxmins=sort([Maxs Mins]);
%     disp('Code to be completed')
%     if isempty(maxmins)
%         maxmins = 1;
%     else
%         if maxmins(1)~=1 && maxmins(end)~=N
%             if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
%                 maxmins=[1 maxmins];
%             end
%         end
%     end

if sum(Maxs==0)>0
    Maxs(Maxs==0)=1;
end
if sum(Mins==0)>0
    Mins(Mins==0)=1;
end
if sum(maxmins==0)>0
    maxmins(maxmins==0)=1;
end

if nargout<=1
    varargout{1}=maxmins;
elseif nargout==2
    varargout{1}=Maxs;
    varargout{2}=Mins;
end

end