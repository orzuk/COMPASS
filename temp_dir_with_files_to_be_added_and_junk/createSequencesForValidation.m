if 1==2 % prepare the files
  clear

  load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human


  % for each amplified sequence - find all sequences in 750 and 450 and 1500 that match it  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1500
  clear all
  global Header_750 Sequence_750 header_uni1to350_primers450  seq_uni1to350_primers450 Header_uni Sequence_uni
  load ~/CS/BAC/primers750/bac16s_primers750_full_without_ambiguous
  load  ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450.mat
  load ~/CS/BAC/bacteria_s16_data_uni

  clear Sequence_amp* Header_amp*

  load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human

  primerName = '205914'

  currSeq = ['Sequence_amp_',primerName];
  currSeq = eval(currSeq);

  currHeader = ['Header_amp_',primerName];
  currHeader = eval(currHeader);

  clear res
  for i=1:length(currSeq)
    res{i}{1}.amp.Sequence = currSeq{i};
    res{i}{1}.amp.Header = currHeader{i};
    for j=1:length(currHeader{i})
      [i,j]
      res{i}{j}.s750 = struct;
      res{i}{j}.s450 = struct;
      [res{i}{j}.s750.Header,res{i}{j}.s750.Sequence,res{i}{j}.s450.Header,res{i}{j}.s450.Sequence,res{i}{j}.Full.Header,res{i}{j}.Full.Sequence]=findPositionOfAmpSequence(currSeq{i},currHeader{i}{j});
    end
  end

  % align and add numbers
  for i=1:length(currSeq)
    
    
    
    for j=1:length(res{i})
      
      currNumber_Full = res{i}{j}.Full.Header;
      a = find(currNumber_Full==' ');
      currNumber_Full = currNumber_Full(1:a(1)-1);
      currNumber_Full = str2num(currNumber_Full);
      
      seq1 = res{i}{1}.amp.Sequence;
      
      % Full
      seq2 = res{i}{j}.Full.Sequence;
      AlignStruct = localalign(seq1, seq2,'Alphabet','NT');
      res{i}{j}.Full.Start = AlignStruct.Start;
      res{i}{j}.Full.Stop = AlignStruct.Stop;
      res{i}{j}.Full.Readme = 'the first is the amplified sequence';
      res{i}{j}.Full.number = currNumber_Full;
      
      
      % align 750 and 450 to full
      seq1Full = res{i}{j}.Full.Sequence;
      
      % s750
      if ~isempty(res{i}{j}.s750.Header)
        seq2 = res{i}{j}.s750.Sequence;
        AlignStruct = localalign(seq1Full, seq2,'Alphabet','NT');
        res{i}{j}.s750.Start = AlignStruct.Start;
        res{i}{j}.s750.Stop = AlignStruct.Stop;
        res{i}{j}.s750.Readme = 'the first is the full sequence';
        currNumber_750 = res{i}{j}.s750.Header;
        a = find(currNumber_750==' ');
        currNumber_750 = currNumber_750(1:a(1)-1);
        currNumber_750 = str2num(currNumber_750);
        res{i}{j}.s750.number = [currNumber_750,currNumber_Full];
        res{i}{j}.s750.numberReadme = 'the first is number in 450 the second in number in full';
        
        
      end % end 750
      
      % 450
      if ~isempty(res{i}{j}.s450.Header)
        seq2 = res{i}{j}.s450.Sequence;
        AlignStruct = localalign(seq1Full, seq2,'Alphabet','NT');
        res{i}{j}.s450.Start = AlignStruct.Start;
        res{i}{j}.s450.Stop = AlignStruct.Stop;
        res{i}{j}.s450.Readme = 'the first is the full sequence';
        
        currNumber_450 = res{i}{j}.s450.Header;
        a = find(currNumber_450==' ');
        currNumber_450 = currNumber_450(1:a(1)-1);
        currNumber_450 = str2num(currNumber_450);
        res{i}{j}.s450.number = [currNumber_450,currNumber_Full];
        res{i}{j}.s450.numberReadme = 'the first is number in 450 the second in number in full';
      end %
    end
    
  end


  save(['~/CS/BAC/12Samples/validation/Human_scf_files/res_',primerName],'res') 
  % add numbers

end
%%%%%%%%%%%%% 21

%%%%%%%%%%% 51
if 1==2 % prepare the files
  clear

  load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human


  % for each amplified sequence - find all sequences in 750 and 450 and 1500 that match it  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1500
  clear all
  global Header_750 Sequence_750 header_uni1to350_primers450  seq_uni1to350_primers450 Header_uni Sequence_uni
  load ~/CS/BAC/primers750/bac16s_primers750_full_without_ambiguous
  load  ~/CS/BAC/primers_startAs750_length350/bac16s_1to350_primers450.mat
  load ~/CS/BAC/bacteria_s16_data_uni

  clear Sequence_amp* Header_amp*

  load ~/CS/BAC/12Samples/validation/Human_scf_files/insilocopcr_sanger_human

  primerName = '205914'

  currSeq = ['Sequence_amp_',primerName];
  currSeq = eval(currSeq);

  currHeader = ['Header_amp_',primerName];
  currHeader = eval(currHeader);

  clear res
  for i=1:length(currSeq)
    res{i}{1}.amp.Sequence = currSeq{i};
    res{i}{1}.amp.Header = currHeader{i};
    for j=1:length(currHeader{i})
      [i,j]
      res{i}{j}.s750 = struct;
      res{i}{j}.s450 = struct;
      [res{i}{j}.s750.Header,res{i}{j}.s750.Sequence,res{i}{j}.s450.Header,res{i}{j}.s450.Sequence,res{i}{j}.Full.Header,res{i}{j}.Full.Sequence]=findPositionOfAmpSequence(currSeq{i},currHeader{i}{j});
    end
  end

  % align and add numbers
  for i=1:length(currSeq)
    
    
    
    for j=1:length(res{i})
      
      currNumber_Full = res{i}{j}.Full.Header;
      a = find(currNumber_Full==' ');
      currNumber_Full = currNumber_Full(1:a(1)-1);
      currNumber_Full = str2num(currNumber_Full);
      
      seq1 = res{i}{1}.amp.Sequence;
      
      % Full
      seq2 = res{i}{j}.Full.Sequence;
      AlignStruct = localalign(seq1, seq2,'Alphabet','NT');
      res{i}{j}.Full.Start = AlignStruct.Start;
      res{i}{j}.Full.Stop = AlignStruct.Stop;
      res{i}{j}.Full.Readme = 'the first is the amplified sequence';
      res{i}{j}.Full.number = currNumber_Full;
      
      
      % align 750 and 450 to full
      seq1Full = res{i}{j}.Full.Sequence;
      
      % s750
      if ~isempty(res{i}{j}.s750.Header)
        seq2 = res{i}{j}.s750.Sequence;
        AlignStruct = localalign(seq1Full, seq2,'Alphabet','NT');
        res{i}{j}.s750.Start = AlignStruct.Start;
        res{i}{j}.s750.Stop = AlignStruct.Stop;
        res{i}{j}.s750.Readme = 'the first is the full sequence';
        currNumber_750 = res{i}{j}.s750.Header;
        a = find(currNumber_750==' ');
        currNumber_750 = currNumber_750(1:a(1)-1);
        currNumber_750 = str2num(currNumber_750);
        res{i}{j}.s750.number = [currNumber_750,currNumber_Full];
        res{i}{j}.s750.numberReadme = 'the first is number in 450 the second in number in full';
        
        
      end % end 750
      
      % 450
      if ~isempty(res{i}{j}.s450.Header)
        seq2 = res{i}{j}.s450.Sequence;
        AlignStruct = localalign(seq1Full, seq2,'Alphabet','NT');
        res{i}{j}.s450.Start = AlignStruct.Start;
        res{i}{j}.s450.Stop = AlignStruct.Stop;
        res{i}{j}.s450.Readme = 'the first is the full sequence';
        
        currNumber_450 = res{i}{j}.s450.Header;
        a = find(currNumber_450==' ');
        currNumber_450 = currNumber_450(1:a(1)-1);
        currNumber_450 = str2num(currNumber_450);
        res{i}{j}.s450.number = [currNumber_450,currNumber_Full];
        res{i}{j}.s450.numberReadme = 'the first is number in 450 the second in number in full';
      end %
    end
    
  end


  save(['~/CS/BAC/12Samples/validation/Human_scf_files/res_',primerName],'res') 
  % add numbers

end
%%%%%%%%% end 51



% do 21
fileName = '5103_O721-F.scf';
SangerSeq_start = 30; % start of Sanger major allele
AmplifiedSeq_start = 1;
lengthValidSanger = 750-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 700]
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

fileName = '5103_O721-R.scf';
SangerSeq_start = 40; % start of Sanger major allele
AmplifiedSeq_start = 1;
lengthValidSanger = 800-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 700]
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

fileName = '5103_O1021-F.scf';
SangerSeq_start = 30; % start of Sanger major allele
AmplifiedSeq_start = 1;
lengthValidSanger = 300-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 200]
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)

fileName = '5103_S721-F.scf';
AmplifiedSeq_start = 100;
SangerSeq_start = 30; % start of Sanger major allele
lengthValidSanger = 800-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 400];
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)
% did not work

fileName = '5103_S1021-F.scf';
AmplifiedSeq_start = 100;
SangerSeq_start = 30; % start of Sanger major allele
lengthValidSanger = 400-SangerSeq_start; % 750 is the length of the Sanger
pos_xl = [100 400];
results_chroma2_21(fileName,res,SangerSeq_start,AmplifiedSeq_start,lengthValidSanger,pos_xl)
% did not work


