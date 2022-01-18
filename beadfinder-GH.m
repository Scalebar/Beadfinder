
%assumptions - all images in directory were taken under the same conditions
%and beads are well spaced. Extracts bead boxes and saves as matrices

Directory="/home/ian/Desktop/psf-see/";
Directoryout="/home/ian/Desktop/psf-see/out";
Beadsize=0.1; %microns

Dirlist=dir(fullfile(Directory,'*.czi'));
Numberoffiles=length(Dirlist);

%Read in first file to get initalisation details
Filenamein=fullfile(Dirlist(1).folder,Dirlist(1).name);

%bfopen - bioformats plugin
Newimage=bfopen(char(Filenamein)); % Zeiss confocal Czi format - not merely images

%Extract two important data tables
Filedetails=Newimage{1,1}{1,2};
Globalsettings=Newimage{1,2};

%Extract relevant data
Xscale=str2num(Globalsettings.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'))*1000000;
Yscale=str2num(Globalsettings.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'))*1000000;
Zscale=str2num(Globalsettings.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1'))*1000000;
% Scale in microns per pixel

Objective=Globalsettings.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|Objective #1');
Numberofchannels=str2num(extractAfter(Filedetails,length(Filedetails)-1));%There's got to be a better way
FrameDimensions=size(Newimage{1,1}{1,1});
Xdim=FrameDimensions(1);
Ydim=FrameDimensions(2);

%Bead box size
Zrange=35;
Xrange=25;
Yrange=25; %Half length of axis for bead box - measuring from centre outwards

% Make a model bead for psf estimation
Xbeadsize=Beadsize/Xscale;
Ybeadsize=Beadsize/Yscale;
Zbeadsize=Beadsize/Zscale; %For current settings this is a dot
Fakebead=zeros(1+(Xrange*2),1+(Yrange*2),1+(Zrange*2)); % puts a pixel in the middle of the volume
Fakebead(Xrange+1,Yrange+1,Zrange+1)=1;



%Start file loop
for Filecounter=1:Numberoffiles
    Filenamein=fullfile(Dirlist(Filecounter).folder,Dirlist(Filecounter).name);
    Filebase=Dirlist(Filecounter).name;

    %bfopen - bioformats plugin
    Newimage=bfopen(char(Filenamein)); % Zeiss confocal Czi format
    Stackdepth=size(Newimage{1});
    Stackdepth=(Stackdepth(1));
    Stackdepth=Stackdepth/Numberofchannels;  % get Z size in a roundabout way
    Channelarray=uint16(zeros(Xdim,Ydim,Stackdepth));

    %Loop through channels
    for Currentchannel=1:Numberofchannels

        % load data from czi format file to matrix
        for Planenumber=Currentchannel:Numberofchannels:Stackdepth
            Framenumber=floor(Planenumber/Numberofchannels)+1;
            Channel=Newimage{1, 1}{Planenumber, 1};
            Channelarray(:,:,Framenumber)=Channel;

        end

        %Generate denoise via blur, binarise to get centroids
        Blurchannel=medfilt3(Channelarray,[5 5 5]);
        Blurthresh=max(max(mean(Channelarray)))/2; %tinkered with max and mean to get right balance
        Binchannel=Blurchannel>Blurthresh;

        %Avoid beads close to edge by adding margin
        Binchannel(1:Xrange,:,:)=0;
        Binchannel(Xdim-Xrange:Xdim,:,:)=0;
        Binchannel(:,1:Yrange,:)=0;
        Binchannel(:,Ydim-Yrange:Ydim,:)=0;
        Binchannel(:,:,1:Zrange)=0;
        Binchannel(:,:,(Stackdepth)-Zrange:(Stackdepth))=0;

        %Apply centroid labels
        Beadlist=regionprops(Binchannel,'Centroid');

        %Extract bead images
        Totalbeads=length(Beadlist);
        Beadarray=zeros(1+(Xrange*2),1+(Yrange*2),1+(Zrange*2)); % Temp array for each image
        %Beadsum=zeros(1+(Xrange*2),1+(Yrange*2),1+(Zrange*2)); %Cumulative total for average

        for Beadcounter=1:Totalbeads
            Beadlocation=round(Beadlist(Beadcounter).Centroid);
            Xstart=Beadlocation(1)-Xrange;
            Xfinish=Xrange+Beadlocation(1);
            Ystart=Beadlocation(2)-Yrange;
            Yfinish=Yrange+Beadlocation(2);

            Beadarray=double(Channelarray(Ystart:Yfinish,Xstart:Xfinish,1:Stackdepth));
            Beadcounter % Just to check progress print the counter
            Beadstring=int2str(Beadcounter);
            Channelstring=int2str(Currentchannel);
            Filenameout=erase(Filebase,".czi")+"-"+Channelstring+"-"+Beadstring+".mat";
            Filenameout=fullfile(Directoryout,Filenameout);
            Beadcounter
            save(Filenameout,'Beadarray');
        end


    end
    
end
Filenameout=erase(Filebase,".czi")+"-"+Channelstring+"-modelbead.mat";
    Filenameout=fullfile(Directoryout,Filenameout);

    save(Filenameout,'Fakebead');





    
    
