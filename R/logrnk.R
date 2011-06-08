####performs log rank test
logrnk<-function(dataL, dataS){
    ####add a new row,and merge all necessary fields in data frames and new row into one dataframe
    new.row<-c(0, FALSE, 0, 0)
    columns_L<-data.frame(cbind(dataL$PatientOrderValidation, as.character(dataL$group), dataL$True_STs, dataL$censored), stringsAsFactors=FALSE)
    columns_S<-data.frame(cbind(dataS$PatientOrderValidation, as.character(dataS$group), dataS$True_STs, dataS$censored), stringsAsFactors=FALSE)
    data.new<-rbind(columns_L, columns_S)
    data.new<-rbind(new.row, data.new)
    names(data.new)<-c("patient", "group", "True_STs", "censored")
    data.new$patient<-as.integer(data.new$patient)
    data.new$group<-as.character(data.new$group)
    data.new$True_STs<-as.numeric(data.new$True_STs)
    data.new$censored<-as.integer(data.new$censored)
    no.patients<-nrow(data.new)-1
    #####order the data in ascending order of survival times, then censorship status
    data.new<-data.new[order((as.numeric(data.new$True_STs)), data.new$censored), ]
    total.risk<-c(no.patients, no.patients:1) #vector of all patients at risk from survivaltimes 1:end
    risk_L<-c(nrow(dataL)) #risk of patients in dataset dataL
    deaths<-c(0) #vector of number of deaths occuring at each time point
    expected_L<-c(0) #vector of the expected number of deaths occuring in group L at each time point
    ####calculate the number of patients in group L who are at risk at each time point
    for(disk in 2:(no.patients + 1)) {
        if(data.new$group[disk-1] == "L") 
            risk_L<-c(risk_L,(risk_L[disk-1] - 1))
        else 
            risk_L<-c(risk_L, risk_L[disk-1])
        if(data.new$censored[disk] == 1)
            mi<-0
        else 
            mi<-1
    deaths<-c(deaths, mi)
    }
    total_deaths<-sum(deaths) ####total number of death events in both data sets
    #####calculate the expected fraction of patients at risk in group L from the total number of patients at risk
    for(patient in 2:length(data.new$group)) {
        expected_L<-c(expected_L, (risk_L[patient]/total.risk[patient]) * deaths[patient])
    }
expected_L_fin<-c(expected_L[length(data.new$group)])
    expected_L_fin<-expected_L[(no.patients + 1)]
    ####sum up the expected values in group L for the same time point
    for(riskL in c(no.patients:1)) {
        if(data.new$True_STs[riskL] != data.new$True_STs[riskL+1])
            expected_L_fin<-c(expected_L_fin, expected_L[riskL])
        else
            expected_L_fin<-c(expected_L_fin, (expected_L_fin[length(expected_L_fin)] + expected_L[riskL]))
    }
    expected_L_fin<-expected_L_fin[length(expected_L_fin):1]
    total_expected_L<-0 ####total number of patients at risk in group L summed over all time points
    deaths_L<-0 #####total number of observed deaths in group L
    for (death in 2:length(data.new$group)) {
        deaths_i<-sum(deaths[death:length(deaths)])
        if (data.new$group[death] == "L") 
            deaths_L<-deaths_L + deaths[death] 
        if(0<deaths_i) {
            if(data.new$True_STs[death-1] != data.new$True_STs[death])
                total_expected_L<-(total_expected_L + expected_L_fin[death])
        }
        total_expected_S<-total_deaths - total_expected_L
        deaths_S<-total_deaths - deaths_L
        Xsq<-(((deaths_L - total_expected_L)^2/total_expected_L) + ((deaths_S - total_expected_S)^2/total_expected_S))
    }
    p_value<-1 - pchisq(Xsq, df=1) ####df=the number of groups-1
    return(list(Xsq=Xsq, pValue=p_value))
}