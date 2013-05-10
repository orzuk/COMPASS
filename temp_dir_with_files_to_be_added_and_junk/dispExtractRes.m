function dispExtractRes(auxData,userDir,nb)

basicDir = [userDir,'CS/BAC/',auxData.outFileDirName,'/'];
for readLength=auxData.readLength
  for np=auxData.pwr
    for nr=auxData.nr
      for addNoiseFlagInput=auxData.addNoiseFlagInput
        for correctMeasurementForReadErrorsFlagInput=auxData.correctMeasurementForReadErrorsFlagInput
          if ~(addNoiseFlagInput==0 & correctMeasurementForReadErrorsFlagInput==1)
            res(readLength,np*10+1,:) = extractRes(basicDir,auxData.outFileDirName,0,nb,np,nr,addNoiseFlagInput,correctMeasurementForReadErrorsFlagInput,readLength,410849,1);
          end
        end
      end
    end
  end
end
