T=readtable('Anti-DE datasets.csv','ReadVariableNames',false);
T.Properties.VariableNames = {'sample_id', 'species','tissue','geo_acc','geo_mat','batch_id'};

% Access list of GEO matrices
acclist = string(T.geo_mat);
sourcelist = string(T.geo_acc);
samplelist = string(T.batch_id);
specieslist = string(T.species);

fid = fopen('sample_cellnum.txt','w');
for k=1:length(samplelist)

    try
        sampleid=samplelist(k);
        sourceid=sourcelist(k);
        specieid=specieslist(k);
        fprintf(" Sample %d with ID %s from %s \n",k,sampleid, sourceid);
        fprintf(" BatchID %s for species %s \n",sampleid, specieid);
        mkdir(sourceid);
        outfile = sprintf('./%s/%s.mat',sourceid,sampleid);
        if exist(outfile,"file")
            load(outfile);
            fprintf(fid,'%d\t%s\t%d\n',k, sampleid,sce.NumCells);
            continue;
        end
        acc = acclist(k);
        [sce] = sc_readgeoaccession(acc);
        metainfo = sprintf("Source: %s - %s", sourceid,sampleid);
        sce = sce.appendmetainfo(metainfo);
        sce.c_batch_id = repmat(sampleid, sce.NumCells, 1);
        sce = sce.qcfilter;
        if sce.NumCells > 0
            sce = sce.embedcells('umap', true, false, 3);
            sce = sce.clustercells([], [], true);
            sce = sce.assigncelltype(specieid, false);
            [c,cL] = grp2idx(sce.c_cell_type_tx);
            sce.c = c; 
        end
        save(outfile,'sce','-v7.3');
    catch ME
        disp(ME.message)
    end
    fprintf(fid,'%d\t%s\t%d\n',k, sampleid,sce.NumCells);
    pause(3);
    fprintf(" \n\n");

end
fclose(fid);
