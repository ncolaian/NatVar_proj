#! /usr/bin/env perl

# This script will run clarks scripts on the cluster and perform the analysis necessary

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;

#My Variables
my $help = 0;
my $man = 0;
my $sample_names; #make sure the names are in the order of your sample
my $out_dir;
my $ros_out;
my $ros_name;
my $to_r = 0;
my $compare_pep_col;
#need to change this if using other machines
my $create_auc_graphs_and_comp_path = "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/";

#Read in the variables from the command line
GetOptions( 'man'   =>  \$man,
            'help'  =>  \$help,
            'sample_names|s=s' =>   \$sample_names,
            'out_dir|o=s'   =>  \$out_dir,
            'ros_data|r=s'  =>  \$ros_out,
            'ros_name|n=s'  =>  \$ros_name,
            'to_r|tr'       =>  \$to_r,
            'pep_num|pn=i'  =>  \$compare_pep_col,
            ) || die("There was an error in the command line arguements\n");

# Pod Usage for the manual and the help pages
if ( $help ) { pod2usage(0) }
if ( $man ) { pod2usage(-verbose =>3) }

# Setup logging environment
my $logger = get_logger();

## Main ##
my $names_aref = format_data();
if ( $to_r == 0) {
    `cd $out_dir`;
    #create_and_run_DataParseAnalyzeCP($names_aref);
    #create_and_run_ROScomparisons($names_aref);
    get_curve_parameters($names_aref);
}
else {
    run_r_script($names_aref);
}


## Subroutines ##
sub format_data {
    if ( !(-d $out_dir) ) {
        mkdir $out_dir;
    }
    open my $IN, "<", $ros_out;
    open my $NAME, "<", $sample_names;
    
    my @names;
    my %data;
    while ( <$NAME> ) {
        chomp $_;
        if ( $_ eq "" ) {
            next;
        }
        push @names, $_;
        $data{$_} = [];
    }
    
    #read in data to hash
    my @time_points;
    my @lines;
    while ( <$IN> ) {
        chomp $_;
        next if ( $_ =~ /^#/ || $_ =~ /^Plate/ || $_ =~ /^Time/ || $_ eq "" || $_ =~ /^,,,,/ || $_ =~ /^$/);
        last if ( $_ =~ /End/);
        #check for lines with times
        #get rid of carraige return
        $_ =~ s/\r//g;
        my @split_line = split /\t/, $_;
        next if ( scalar @split_line == 0);
        if ( $split_line[0] =~ /:/ ) {
			if (scalar(@lines) > 0) {
				for ( my $i = 0; $i< scalar(@names); $i++) {
					push $data{$names[$i]}, $lines[$i];
				}
			}
            @lines = ();
            $split_line[0] = (split /:/, $split_line[0])[0] * 60.000;
            push @time_points, shift @split_line;
            shift @split_line; #throw away temp
            push @lines, \@split_line;
            #go through and save the data for that time point
        }
        else {
            shift @split_line;
            shift @split_line;
            push @lines, \@split_line;
        }
        
    }
    #get last line
    for ( my $i = 0; $i< scalar(@names); $i++) {
        push $data{$names[$i]}, $lines[$i];
    }
    
    close $NAME;
    close $IN;
    #print out data
    open ( my $OUT, ">", "$out_dir/$ros_name.txt" );
    #foreach my $name ( @names ) {
    #    print $OUT "$name\nTime0;A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;\n";
    #    my $mm_name_aref = $data{$name};
    #    for ( my $i = 0; $i < scalar(@time_points); $i++ ) {
    #        my $time = $time_points[$i];
    #        print $OUT "$time;" . join(";", @{$mm_name_aref->[$i]}) . ";\n";
    #    }
    #    print $OUT "\n";
    #}
    #else {
        print $OUT "Name,Time,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12\n";
        foreach my $name ( @names ) {
            my $mm_name_aref = $data{$name};
            for ( my $i = 0; $i < scalar(@time_points); $i++ ) {
                my $time = $time_points[$i];
                print $OUT "$name,$time," . join(",", @{$mm_name_aref->[$i]}) . "\n";
            }
        }
    close $OUT;
    return(\@names);
}

sub run_r_script {
    my ($names_aref) = @_;
    my $compare_pep = $names_aref->[$compare_pep_col-1];
    my $cmd = "Rscript $create_auc_graphs_and_comp_path/create_auc_graphs_and_comp.R -o $out_dir -b $ros_name -i $out_dir/$ros_name.txt -r $compare_pep";
    `$cmd`;
    return 1;
}

sub get_curve_parameters {
    my ( $names_aref ) = @_;
    my $cmd = "Rscript $create_auc_graphs_and_comp_path/calculate_curve_params.R -o $out_dir -i $out_dir/$ros_name.txt -n $ros_name";
    my $sbatch = "sbatch -t 24:00:00 -o $out_dir/curve_params.out --open-mode=append --mem=8000 -J $ros_name.curve --wrap=\"$cmd\"";
    system($sbatch);
    return(1);
}

sub sumbmit_and_stall {
    my ( $job_aref, $job_name, $out_directory, $mem, $time ) = @_;
    
    foreach my $cmd ( @$job_aref ) {
        my $sbatch_cmd = "sbatch -t $time -o $out_directory/$job_name.out --open-mode=append --mem=$mem -J $job_name --wrap=\"$cmd\"";
        system($sbatch_cmd);
    }
    
    #Create Active job file
    `echo start > $out_directory/ACTIVE_JOBS`;
    
    my $minutes = 0;
    while ( -s "$out_directory/ACTIVE_JOBS" ) {
        sleep 60;
        `squeue --name=$job_name -h > $out_directory/ACTIVE_JOBS`;
        $minutes++;
    }
    
    `rm $out_directory/ACTIVE_JOBS`;
    $logger->info("All $job_name jobs have finished in ~$minutes minutes");
    
    return 1;
}

#### DEPRICATED ####

sub create_and_run_DataParseAnalyzeCP {
    my ($names_aref) = @_;
    open my $OUT, ">", "$out_dir/DPACP.m";
    my $file = "cd $out_dir; \n
fid = fopen('$ros_name');
m=[];
numExperiments=0;
while 1
tline = fgetl(fid);
if ~ischar(tline)
break
end

tline;

saveFiles=1;

if(saveFiles==1)
fileNames={";

foreach my $name ( @$names_aref ) {
    $file .= "'$name',";
}
chop $file;
$file .= "}
end;

s=regexpi(tline, '[0-9_]+[0-9_]+:+[0-9_]+[0-9]+:+[0-9_]+[0-9]','match');
[mat tok ext] = regexp(tline, '[0-9]*[,.]?[0-9]+;', 'match', 'tokens', 'tokenExtents');
%mat
%pause
d=[];
if(~isempty(mat))
    numTimePoints=size(mat,2);
    for(i=1:size(mat,2))
        i
        class(mat(i))
        [mat2 tok ext] = regexp(cell2mat(mat(i)), '[0-9]*[,.]?[0-9]+;', 'match','tokens', 'tokenExtents');
        mat2
        d= [d, str2num(cell2mat(mat2))];
    end
        
end
m=[m;d]
end
m(:,1)=[]; % first column is zeros.... rip these off


%m currently contains all experiments.
%parse out each experiment and put into an indexed matrix
mm=[];
counter=0;
for(i=1:size(m,1))
    %m(i,:)
    %sum(m(i,:))
    %pause
    if(sum(m(i,:))==sum(1:numTimePoints-1))
        fprintf('new experiment seen\\n')
        mm
        
        if(numExperiments>0)
            experimentsHold(:,:,numExperiments) = mm;
        end
        
        mm=[];
        numExperiments=numExperiments+1;
        fprintf('-------------------\\n') 

    else
        mm=[mm;m(i,:)]; 
    end
end
experimentsHold(:,:,numExperiments) = mm;

[timePts,reps]=size(experimentsHold(:,:,1))

%-------- Plot the experiments-----------------
figure(1)
clf
hold on
%for(i=1:numExperiments)
for(i=1:1)
i
    %pause
    col=[rand,rand,rand];
    %for(j=1:timePts)
    for(k=1:reps)
        k
        %pause
        col=[rand,rand,rand];
        %rep=experimentsHold(j,:,i)
        rep=experimentsHold(:,k,i)
        i*ones(size(rep))
        plot(1:timePts,rep','o-','COLOR',col)
        %plot(j*ones(size(rep)),rep,'.','COLOR',col)
    end
end
%pause(1)

%----------------------------- Everything above is just a parser and data
%%-----------------------------formatter.
%%%-----------------------------Below is the curve fitter.
FID = fopen('$out_dir/$ros_name.DPAout.txt','w')
if(saveFiles==1)
if(length(fileNames)~=numExperiments)
   error('filenames inconsistent with the number of experiments\\n') 
end
end

for(iii=1:numExperiments)

%%fit curve to data
experiment=experimentsHold(:,:,iii);
maxHeight = max(max(experiment));
experiment=experiment/maxHeight; %rescale the experiments so the luminoscity isn't so large
[timePts,reps]=size(experiment)

%form likelihood function

offset=-1;

f2Const = @(theta) theta(1)*exp(-1.*exp(theta(3).*(theta(4)+theta(2)))).*exp(theta(5)*(theta(4)));
logCurve = @(theta,time) theta(1).*exp(-1.*exp(theta(3).*(time+theta(2)))).*((time)<theta(4))+((time)>=theta(4)).*f2Const(theta).*exp(-theta(5)*(time));

timePt=1;
exper=3;
timeScale=repmat(1:timePts,reps,1)';

[maxE,ind]=max(experiment); %set prior bounds for change point. % this really only accelerates convergence
maxTime=max(ind)+2;
minTime=min(ind)-2;

%minTime=3;
%maxTime=6;


logpost = @(params)log(params(6))*(timePts*reps-2)/2-params(6)/2*sum(sum( (experiment-logCurve([params(1),params(2),params(3),params(4),params(5)],timeScale)).^2 ))...
+ log(double(params(1)>0))+ log(double(params(1)<3.2))+log(double(params(3)<-2))+log(double(params(4)>minTime))+log(double(params(4)<maxTime))+log(double(params(5)>0))... %+log(double(params(2)<0))+log(double(params(3)<0))+log(double(params(4)>0));
+ log(double(params(3)>-50))+ log(double(params(6)>0))+ log(double(params(2)<0));
%Note:params 1,2,3...5 are theta in logCurve.... theta 6 is phi=1/sigma^2;
logpostH0 = @(params)log(params(6))*(timePts*reps-2)/2-params(6)/2*sum(sum( (experiment-logCurve([params(1),params(2),params(3),params(4),params(5)],timeScale)).^2 ))...
+ log(double(params(1)>0))+ log(double(params(1)<3.2))+log(double(params(3)<0))+log(double(params(3)>-.001))+log(double(params(4)>minTime))+log(double(params(4)<maxTime))+log(double(params(5)>0))+log(double(params(5)<0.001))... %+log(double(params(2)<0))+log(double(params(3)<0))+log(double(params(4)>0));
+ log(double(params(6)>0))+ log(double(params(2)<0))+ log(double(params(2)>-.05));

initial = [.1, -.1, -4,maxTime/2+minTime/2,1,100.1];
nsamples = 10000;
trace = slicesample(initial,nsamples,'logpdf',logpost,'width',[5, 5, 20,5, 5,20],'thin',20,'burnin',1000);
initialH0 = [.1, -.001, -.0001,maxTime/2+minTime/2,.0001,.1];
traceH0 = slicesample(initialH0,nsamples,'logpdf',logpostH0,'width',[5, .01, .0001,1, .0001,10],'thin',20,'burnin',1000);

traceHold=trace;
traceHold(:,1) = maxHeight.*traceHold(:,1)

if(saveFiles==1)
    fileName=[char(fileNames(iii)),'.mat']
    save(fileName,'traceHold')
end

%plot
figure(2)
clf
k=size(trace,2); %number of parameters.
ii=1;
for(i=1:k)
    subplot(k,2,ii)
    plot(trace(:,i),'r')
    ii=ii+2;
end

k=size(trace,2); %number of parameters.
ii=2
for(i=1:k)
    subplot(k,2,ii)
    plot(traceH0(:,i))
    ii=ii+2;
end

%plot curve over data
h=figure(3)
clf
hold on
%for(i=1:numExperiments)
for(i=1:numExperiments)
i
    %pause
    col=[rand,rand,rand];
    %for(j=1:timePts)
    for(k=1:reps)
        %k
        %pause
        col=[50/255,50/255,50/255]%[rand,rand,rand];
        %rep=experimentsHold(j,:,i)
        rep=experiment(:,k)
        i*ones(size(rep))
        plot(1:timePts,rep','o-','COLOR',col)
        %plot(j*ones(size(rep)),rep,'.','COLOR',col)
    end
end



mt=mean(trace);
ts=1:.1:timePts;
%plot(ts,logCurve(mt(1:6),ts),'k','linewidth',2)

postCurves=zeros(length(ts),nsamples);
for(i=1:nsamples)
  postCurves(:,i)=logCurve(trace(i,:),ts);
  postCurvesH0(:,i)=logCurve(traceH0(i,:),ts);
end
plot(ts,mean(postCurves'),'r','linewidth',2);
plot(ts,mean(postCurvesH0'),'b','linewidth',2);

if(saveFiles==1)
    fileName=[char(fileNames(iii)),'.fig']
    saveas(h,fileName)  %save the figure with the appropriate file name
end

sc=sort(postCurves')';
upperCurve = sc(:,floor(.975*nsamples));
lowerCurve = sc(:,floor(.025*nsamples));
%plot(ts,upperCurve,'b','linewidth',2)
%plot(ts,lowerCurve,'b','linewidth',2)
%get upper and lower 95 for curves

%Compute BIC
med = median(trace);
medH0=median(traceH0);

postHAlt=zeros(length(trace),1);
for(i=1:length(trace))
    postHAlt(i)=logpost(trace(i,:));
end
postH0=zeros(length(traceH0),1);
for(i=1:length(traceH0))
    postH0(i)=logpostH0(traceH0(i,:));
end
%-2log(P(x|theta,k))\approx -2*log(L(B|X))+p*log(N);

%BICNull = -2*logpost(medH0)+4*log(timePts*reps)  %% check if the signs are right.
%BICAlt  = -2*logpost(med)+6*log(timePts*reps)
BICNull = -2*max(postH0)+3*log(timePts*reps)  %% check if the signs are right.
BICAlt  = -2*max(postHAlt)+6*log(timePts*reps)

if(saveFiles==1)
    fileNames(iii)
end

BIC = BICNull-BICAlt;
BF = exp(-0.5*BIC)
post = inv(1+inv(BF))
fprintf(FID,'%d %30s %6.4f   %6.4f\\n',iii,char(fileNames(iii)),BF, post);
for(i=1:1)
    iii
    %pause(1)
    beep
end
%pause
end
fclose(FID)";
    print $OUT $file;
    close $OUT;
    my $cmd = "matlab -r 'DPACP'";
    my @jobs;
    push @jobs, $cmd;
    sumbmit_and_stall(\@jobs, "DPACP", $out_dir, "1000", "24:00:00");
    return 1;
}

sub create_and_run_ROScomparisons {
    my ($name_aref) = @_;
    
    #here we assume that you want to compare everything to the regular and that will be the first value and the negative which is the second value
    my $reg_pep = $names_aref->[$compare_pep_col-1];
    #open an out file
    open my $OUT, ">", "$out_dir/ros_compare.m";
    print $OUT "fid = fopen('$out_dir/$ros_name.comparisons.txt','w');\n";
    
    for (my $i = 0; $i < scalar(@$name_aref); $i++) {
        next if ( $i == $compare_pep_col-1);
        my $pep_name = $name_aref->[$i];
        my $out_line = "data1=load('$reg_pep.mat');
data2 = load('$pep_name.mat');

tr1=data1.traceHold;
tr2=data2.traceHold;

traceDiff=tr1-tr2;

%plot
figure(2)
clf
k=size(tr1,2); %number of parameters.
ii=1;
for(i=1:k)
    subplot(k,1,i)
    plot(traceDiff(:,i),'r')
end

for(i=1:6)
    m(i)=mean(traceDiff(:,i));
    sortTrace = sort(traceDiff(:,i));
    lowerLimit(i) = sortTrace(ceil(0.05*length(sortTrace)));
    upperLimit(i) = sortTrace(floor(0.95*length(sortTrace)));
end



fprintf(fid,'$reg_pep       , $pep_name       , ')
for(i=1:5)
    fprintf(fid,'%6.2f ,',m(i))
end
fprintf(fid,'\\n')
fprintf(fid,'                          ,                        , ')
for(i=1:5)
    fprintf(fid,'[%6.2f, %6.2f] ,',lowerLimit(i),upperLimit(i))
end
fprintf(fid,'\\n')
";
        print $OUT $out_line;
        
    }
    print $OUT "\nfclose(fid)";
    close $OUT;
    my $cmd = "matlab -r 'ros_compare'";
    my @jobs;
    push @jobs, $cmd;
    sumbmit_and_stall(\@jobs, "compare", $out_dir, "1000", "24:00:00");
    return 1;
}

__END__
=head1 TITLE



=head1 VERSION



=head1 INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Readonly;
Path::Class;
Data::Dumper;
Log::Log4perl qw(:easy);
Log::Log4perl::CommandLine qw(:all);

=head1 INHERIT

=head1 SYNOPSIS

=head1 PARAMETERS

=head1 CONFIGURATION AND ENVIRONMENT

    Need to load muscle and BLAST

=head1 DEPENDENCIES

    muscle
    BLAST
    Master_aln
    
=head1 INCOMPATIBILITIES

    None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests	
	
=head1 AUTHOR

Nicholas Colaianni
contact via C<< <ncolaian@live.unc.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2017, Nicholas Colaianni
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut
