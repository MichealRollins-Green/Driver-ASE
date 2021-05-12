function OutPutDriverParametersIntoTxt(dsd,outfilepath)
d=dsd;
out.data=[d.p d.fdr d.fm];
out.rowlabels=d.rowlabels;
out.collabels={'-log10P','fdr','fm'};
print_excel_sheet(out,outfilepath);