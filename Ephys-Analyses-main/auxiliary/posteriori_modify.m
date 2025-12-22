for subj=1:n_subjects
    path = strcat(directory_name,os_flag,subj_list{subj});
    results = load(path);
    %results.LFP_results.results.coherency.baseline.PFC_rs.freq = 0:0.1:100
    %results.LFP_results.results.coherency.cno.PFC_rs.freq = 0:0.1:100
    if isfield(results.MUA_results.time,'CNO')
        results.MUA_results.time.cno = results.MUA_results.time.CNO;
        results.MUA_results.time = rmfield(results.MUA_results.time,'CNO')
        MUA_results = results.MUA_results
        %save(path,'LFP_results','-v7.3')
        save(path,'MUA_results','-v7.3')
    end

end