import java.io.*;
import java.sql.SQLOutput;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ACAT {

    public static int count;
    public static String Gene_ENSG_ID;
    public static int Gene_ENSG_Start;
    public static int Gene_ENSG_End;
    public static String Gene_Chromosome;

    public static HashMap<String, String> SNP_EN_Weights_Dic;
    public static HashMap<String, String> EN_SNPs_RSID_Dic;
    public static HashMap<String, String> EN_SNPs_Ref_Dic;
    public static HashMap<String, String> EN_SNPs_Alt_Dic;

    public ACAT(){}

    public static void get_en_info(String weight_file) throws IOException {
        BufferedReader bufferedReader=new BufferedReader(new FileReader(weight_file));
        String line=bufferedReader.readLine();
        boolean flag=false;
        while(line!=null){
            String[] tmp;
            line = line.trim();
            tmp=line.split("\t");
            if (Gene_ENSG_ID.contains("ENSG")?(tmp[0].split("\\.")[0].equals(Gene_ENSG_ID.split("\\.")[0])):(tmp[0].equals(Gene_ENSG_ID))){
                flag=true;

                String chrome=tmp[1].replace("chr","");
                Gene_Chromosome = chrome;
                Gene_ENSG_Start=Integer.parseInt(tmp[2]);
                Gene_ENSG_End=Integer.parseInt(tmp[3]);
                for(int i = 4; i < tmp.length; ++i) {
                    String[] tmp1 = tmp[i].split("_");
                    try {
                        Double.parseDouble(tmp1[4]);
                    }catch (Exception e){
                        continue;
                    }
                    EN_SNPs_RSID_Dic.put(chrome + ":" + tmp1[1], tmp1[0]);
                    SNP_EN_Weights_Dic.put(chrome + ":" + tmp1[1], tmp1[4]);
                    EN_SNPs_Ref_Dic.put(chrome + ":" + tmp1[1], tmp1[2]);
                    EN_SNPs_Alt_Dic.put(chrome + ":" + tmp1[1], tmp1[3]);
                }
                break;
            }
            line= bufferedReader.readLine();
        }
        bufferedReader.close();

        if(!flag){
            System.out.println("Can not find the gene:\t" + Gene_ENSG_ID + " in the " + weight_file + "!");
        }


    }
    public static String csv2tfile(String genotypeFile, String phenotypeFile, String columns,String output_folder) throws IOException {
        String[] info=columns.split(",");
        int sample_index=Integer.parseInt(info[0]);
        int pheno_index=Integer.parseInt(info[1]);
        BufferedReader br_geno=new BufferedReader(new FileReader(genotypeFile));
        BufferedWriter bw_geno=new BufferedWriter(new FileWriter(output_folder+"plink.tped"));
        String line=br_geno.readLine();
        while(line!=null){
            String[] data=line.trim().split(",");
            String chr=data[0].replace("chr","");
            if(chr.equals(Gene_Chromosome)){
                String index=chr+":"+data[1];
                if(EN_SNPs_Ref_Dic.containsKey(index)){
                    String record=chr+"\t"+EN_SNPs_RSID_Dic.get(index)+"\t0\t"+data[1];
                    String ref=EN_SNPs_Ref_Dic.get(index);
                    String alt=EN_SNPs_Alt_Dic.get(index);
                    for(int i=2;i<data.length;i++){
                        switch (data[i]){
                            case "0":
                                record=record+"\t"+ref+"\t"+ref;
                                break;
                            case "1":
                                record=record+"\t"+ref+"\t"+alt;
                                break;
                            case "2":
                                record=record+"\t"+alt+"\t"+alt;
                                break;
                        }
                    }
                    bw_geno.write(record+"\n");
                }
            }
            line=br_geno.readLine();
        }
        bw_geno.close();
        br_geno.close();

        BufferedReader br_pheno=new BufferedReader(new FileReader(phenotypeFile));
        BufferedWriter bw_pheno=new BufferedWriter(new FileWriter(output_folder+"plink.tfam"));
        br_pheno.readLine();
        String line1=br_pheno.readLine();
        while(line1!=null){
            String[] data=line1.trim().split("\t");
            bw_pheno.write(data[sample_index-1]+"\t"+data[sample_index-1]+"\t0\t0\t0\t"+data[pheno_index-1]+"\n");
            line1=br_pheno.readLine();
        }
        bw_pheno.close();
        br_pheno.close();
        return output_folder+"plink";
    }
    public static String csv2tfile_all(String weightFile,String genotypeFile, String phenotypeFile, String columns,String output_folder) throws IOException {
        //get infomation of snps
        HashMap<String,String> SNPs_Dic=new HashMap<>();
        BufferedReader bufferedReader=new BufferedReader(new FileReader(weightFile));
        String br_line=bufferedReader.readLine();
        while(br_line!=null){
            String[] tmp;
            br_line = br_line.trim();
            tmp=br_line.split("\t");
            String chrome=tmp[1].replace("chr","");
            for(int i = 4; i < tmp.length; ++i) {
                String[] tmp1 = tmp[i].split("_");
                try {
                    Double.parseDouble(tmp1[4]);
                }catch (Exception e){
                    continue;
                }
                SNPs_Dic.put(chrome + ":" + tmp1[1], tmp1[0]+"_"+tmp1[2]+"_"+tmp1[3]);
            }
            br_line= bufferedReader.readLine();
        }
        bufferedReader.close();

        String[] info=columns.split(",");
        int sample_index=Integer.parseInt(info[0]);
        int pheno_index=Integer.parseInt(info[1]);
        BufferedReader br_geno=new BufferedReader(new FileReader(genotypeFile));
        BufferedWriter bw_geno=new BufferedWriter(new FileWriter(output_folder+"plink.tped"));
        String line=br_geno.readLine();
        while(line!=null){
            String[] data=line.trim().split(",");
            String chr=data[0].replace("chr","");
            String index=chr+":"+data[1];
            if(SNPs_Dic.containsKey(index)){
                String[] snp_info=SNPs_Dic.get(index).split("_");
                String record=chr+"\t"+snp_info[0]+"\t0\t"+data[1];
                String ref=snp_info[1];
                String alt=snp_info[2];
                for(int i=2;i<data.length;i++){
                    switch (data[i]){
                        case "0":
                            record=record+"\t"+ref+"\t"+ref;
                            break;
                        case "1":
                            record=record+"\t"+ref+"\t"+alt;
                            break;
                        case "2":
                            record=record+"\t"+alt+"\t"+alt;
                            break;
                    }
                }
                bw_geno.write(record+"\n");
            }
            line=br_geno.readLine();
        }
        bw_geno.close();
        br_geno.close();

        BufferedReader br_pheno=new BufferedReader(new FileReader(phenotypeFile));
        BufferedWriter bw_pheno=new BufferedWriter(new FileWriter(output_folder+"plink.tfam"));
        br_pheno.readLine();
        String line1=br_pheno.readLine();
        while(line1!=null){
            String[] data=line1.trim().split("\t");
            bw_pheno.write(data[sample_index-1]+"\t"+data[sample_index-1]+"\t0\t0\t0\t"+data[pheno_index-1]+"\n");
            line1=br_pheno.readLine();
        }
        bw_pheno.close();
        br_pheno.close();
        return output_folder+"plink";
    }
    public static String plink2summary(String plink,String plink_file, String filetype,String output_folder) throws IOException, InterruptedException {
        ProcessBuilder CMDLine = new ProcessBuilder(plink, "--"+filetype, plink_file,"--assoc","--allow-no-sex","--out", output_folder+"plink");
        Process CMDProcess = CMDLine.start();
        BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));

        String line;
        while((line = br.readLine()) != null) {
            line = line.replace("\n", "");
            System.out.println(line);
        }

        BufferedReader br_error = new BufferedReader(new InputStreamReader(CMDProcess.getErrorStream()));
        while((line = br_error.readLine()) != null) {
            line = line.replace("\n", "");
            System.out.println(line);
        }
        CMDProcess.waitFor();
        String path=output_folder+"plink.assoc";
        File ff=new File(path);
        if (ff.exists())
            return path;
        else
            return output_folder+"plink.qassoc";
    }
    public static void get_weight_p(String summary_file,int window,String columns,String output_folder) throws IOException {
        String[] info=columns.split(",");
        int chr_index=Integer.parseInt(info[0]);
        int pos_index=Integer.parseInt(info[1]);
        int p_index=Integer.parseInt(info[2]);
        BufferedReader bufferedReader=new BufferedReader(new FileReader(summary_file));
        String line=bufferedReader.readLine();
        String file = output_folder + "weight_p.txt";
        BufferedWriter bw = new BufferedWriter(new FileWriter(file, false));
        bw.write("snp\twgt\tp\n");
        line=bufferedReader.readLine();
        while(line!=null){
            // CHR            SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
            // CHR          SNP         BP    NMISS       BETA         SE         R2        T            P
            String[] tmp;
            line = line.trim();
            tmp=line.split("\\s+");
            String chr=tmp[chr_index-1].replace("chr","");
            if(chr.equals(Gene_Chromosome)){
                int pos=Integer.parseInt(tmp[pos_index-1]);
                if((pos>=Gene_ENSG_Start-window) && (pos<=Gene_ENSG_End+window)){
                    String index=chr+":"+tmp[pos_index-1];
                    if(SNP_EN_Weights_Dic.containsKey(index)){
                        //&&tmp[1].equals(EN_SNPs_RSID_Dic.get(index))
                        try{
                            Double.parseDouble(tmp[p_index-1]);
                        }catch (Exception e){
                            line=bufferedReader.readLine();
                            continue;
                        }

                        bw.write(EN_SNPs_RSID_Dic.get(index)+"\t"+SNP_EN_Weights_Dic.get(index).replace("-","")+"\t"+tmp[p_index-1]+"\n");
                        count++;
                    }
                }
            }
            line=bufferedReader.readLine();
        }
        bw.close();
        bufferedReader.close();
    }
    public static void executeACAT(String output_folder,String rscript) throws IOException, InterruptedException {
        String acat_file=output_folder+"ACAT.R";
        BufferedWriter bw = new BufferedWriter(new FileWriter(acat_file, false));
        bw.write("library(ACAT)\n");
        bw.write("library(bigreadr)\n");
        bw.write("library(data.table)\n");
        bw.write("weight_p <- fread(file=\""+output_folder+"weight_p.txt\",sep=\"\\t\")\n");
        bw.write("pval <- ACAT(Pvals=weight_p$p,weights=weight_p$wgt)\n");
        bw.write("result <- data.table(Gene=\""+Gene_ENSG_ID+"\",Pval=pval)\n");
        bw.write("fwrite2(result,file=\""+output_folder+"ACAT.result\",sep=\"\\t\",row.names=F,col.names=T)\n");
        bw.close();

        ProcessBuilder CMDLine = new ProcessBuilder(rscript, output_folder + "ACAT.R");
        Process CMDProcess = CMDLine.start();
        BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));

        String line;
        while((line = br.readLine()) != null) {
            line = line.replace("\n", "");
            System.out.println(line);
        }

        BufferedReader br_error = new BufferedReader(new InputStreamReader(CMDProcess.getErrorStream()));
        boolean flag=false;
        line = br_error.readLine();
        if(line!=null){
            flag=true;
        }
        while( line != null) {
            line = line.replace("\n", "");
            System.out.println(line);
            line=br_error.readLine();
        }
        CMDProcess.waitFor();

        String path = output_folder + "ACAT.result";
        File file=new File(path);
        if (!file.exists()){
            BufferedWriter bw_error=new BufferedWriter(new FileWriter(path));
            bw_error.write("Gene\tPval\n");
            bw_error.write(Gene_ENSG_ID+"\tNA\n");
            bw_error.close();
        }
    }
    public static void get_result(String output_folder) throws IOException{
        String path = output_folder + "ACAT.result";
        BufferedReader br_info = new BufferedReader(new FileReader(path));
        String line="";
        System.out.println("<<The p-value of ACAT on GENE " + Gene_ENSG_ID +">>");
        while((line = br_info.readLine()) != null) {
            line = line.trim();
            String[] tmp = line.split("\t");
            for (String t:tmp) {
                System.out.print(t+"\t");
            }
            System.out.println();
        }
        br_info.close();
        System.out.println("<<==============================================>>");
    }
    public static String get_result1(String output_folder)throws IOException{
        String path = output_folder + "ACAT.result";
        BufferedReader br_info = new BufferedReader(new FileReader(path));
        br_info.readLine();
        String line= br_info.readLine();
        br_info.close();
        return line;
    }



    public static void main(String[] args) throws IOException, InterruptedException {
        if (args.length==0){
            System.out.println("The help information of ACAT:");
            System.out.println("please use");
            System.out.println("------------");
            System.out.println("library(devtools)");
            System.out.println("devtools::install_github(\"yaowuliu/ACAT\")");
            System.out.println("------------");
            System.out.println("install necessary R package");
            System.out.println("-format summary|csv|tfile|bfile|plink The format of input files.");
            System.out.println("-summary The GWAS summary file is used to execute ACAT.");
            System.out.println("-info_columns When the format of input file is summary, user need to choose columns of chromosome, position and p-value which are split by comma,like 1,2,3.");
            System.out.println("-genotype_file When the format of input file is csv, user need to input the file name of genotype.");
            System.out.println("-phenotype_file When the format of input file is csv, user need to input the file name of phenotype.");
            System.out.println("-plink_file When the format of input file is tfile|bfile|plink(tfile representing tped and tfam, bfile representing bim, bed and fam, plink representing ped and map), user need to input the prefix of file name.");
            System.out.println("-weights The weights are used to execute ACAT.");
            System.out.println("-gene The gene is chosen to associate with phenotype.");
            System.out.println("-genelist The list of genes needed to be analyse.");
            System.out.println("-window The window size will be expend around the gene to include SNPs.");
            System.out.println("-Rscript The absolute file path of R.");
            System.out.println("-plink The absolute file path of plink.");
            System.out.println("-output_folder The absolute file path of output folder.");
            System.exit(0);
        }

        String function = args[0];
        if (args[0].equals("ACAT")) {
            DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
            if (args.length < 6) {
                System.out.println("Please input enough parameters to execute ACAT!\t");
                System.exit(0);
            } else {
                String format=null;
                //summary 数据
                String summary_file = null;
                String columns=null;

                //csv数据
                String genotype_file = null;
                String phenotype_file = null;

                //plink数据
                String plink_file=null;

                String weights_file = null;
                String gene_id = null;
                String genelist=null;
                String window_size=null;
                String rscript = null;
                String plink = null;
                String output_folder = null;


//                format="tfile";
//                plink_file="test/CAD";
//                weights_file="test/weight.txt";
//                gene_id="ENSG00000237491";
//                rscript="rrr";
//                plink="utils/plink";
//                output_folder="test/";

                //get input parameters
                for(int k=1;k< args.length;k++){
                    if(args[k].startsWith("-")){
                        switch (args[k]) {
                            case "-format":
                                format = args[k + 1];
                                break;
                            case "-summary":
                                summary_file = args[k + 1];
                                break;
                            case "-info_columns":
                                columns = args[k + 1];
                                break;
                            case "-genotype_file":
                                genotype_file = args[k + 1];
                                break;
                            case "-phenotype_file":
                                phenotype_file = args[k + 1];
                                break;
                            case "-plink_file":
                                plink_file = args[k + 1];
                                break;
                            case "-weights":
                                weights_file = args[k + 1];
                                break;
                            case "-gene":
                                gene_id = args[k + 1];
                                break;
                            case "-genelist":
                                genelist = args[k + 1];
                                break;
                            case "-window":
                                window_size = args[k + 1];
                                break;
                            case "-Rscript":
                                rscript = args[k + 1];
                                break;
                            case "-plink":
                                plink = args[k + 1];
                                break;
                            case "-output_folder":
                                output_folder = args[k + 1];
                                break;
                        }
                    }
                }
                //check inputs
                if (format==null||format.startsWith("-")){
                    System.out.println("Please input the format!");
                    System.exit(0);
                }
                switch (format){
                    case "summary":
                        if (summary_file==null||summary_file.startsWith("-")){
                            System.out.println("Please input the summary file!");
                            System.exit(0);
                        }
                        if (columns==null||columns.startsWith("-")){
                            System.out.println("Please input the info columns!");
                            System.exit(0);
                        }
                        break;
                    case "csv":
                        if (genotype_file==null||genotype_file.startsWith("-")){
                            System.out.println("Please input the genotype file!");
                            System.exit(0);
                        }
                        if (phenotype_file==null||phenotype_file.startsWith("-")){
                            System.out.println("Please input the phenotype file!");
                            System.exit(0);
                        }
                        if (plink==null||plink.startsWith("-")){
                            System.out.println("Please input the plink!");
                            System.exit(0);
                        }
                        break;
                    default:
                        if (plink_file==null||plink_file.startsWith("-")){
                            System.out.println("Please input the plink file!");
                            System.exit(0);
                        }
                        if (plink==null||plink.startsWith("-")){
                            System.out.println("Please input the plink!");
                            System.exit(0);
                        }
                }


                if (weights_file==null||weights_file.startsWith("-")){
                    weights_file="weight.txt";
                }
                boolean geneFlag=true;
                if (gene_id==null||gene_id.startsWith("-")){
                    if (genelist==null||genelist.startsWith("-")){
                        System.out.println("Please input the gene ID!");
                        System.exit(0);
                    }else {
                        geneFlag=false;
                    }
                }
                if (window_size==null||window_size.startsWith("-")){
                    window_size="0";
                }
                if (rscript==null||rscript.startsWith("-")){
                    System.out.println("Please input the Rscript!");
                    System.exit(0);
                }
                if (output_folder==null||output_folder.startsWith("-")){
                    System.out.println("Please input the output folder!");
                    System.exit(0);
                }



                File folder=new File(output_folder);
                if(!folder.exists()){
                    if(!folder.mkdir()){
                        System.out.println("output folder making failed!");
                        System.exit(0);
                    }
                }
                output_folder+="/";



                System.out.println("===============START===============");

                if(geneFlag){
                    count=0;
                    Gene_ENSG_ID=gene_id;
                    EN_SNPs_RSID_Dic=new HashMap<>();
                    SNP_EN_Weights_Dic=new HashMap<>();
                    EN_SNPs_Ref_Dic=new HashMap<>();
                    EN_SNPs_Alt_Dic=new HashMap<>();
                    get_en_info(weights_file);
                    System.out.println("Start to format files.\t"+dtf.format(LocalDateTime.now()));
                    switch (format){
                        case "csv":
                            plink_file=csv2tfile(genotype_file,phenotype_file,"1,2",output_folder);
                            summary_file=plink2summary(plink,plink_file,"tfile",output_folder);
                            columns="1,3,9";
                            break;
                        case "tfile":
                            summary_file=plink2summary(plink,plink_file,"tfile",output_folder);
                            columns="1,3,9";
                            break;
                        case "bfile":
                            summary_file=plink2summary(plink,plink_file,"bfile",output_folder);
                            columns="1,3,9";
                            break;
                        case "plink":
                            summary_file=plink2summary(plink,plink_file,"file",output_folder);
                            columns="1,3,9";
                            break;
                    }
                    get_weight_p(summary_file,Integer.parseInt(window_size),columns,output_folder);
                    System.out.println("End to format files.\t"+dtf.format(LocalDateTime.now()));
                    System.out.println("-------------------------------");
                    if(count==0){
                        System.out.println("There exists no snp around this gene!");
                    }else{
                        System.out.println("Running ACAT with the weights of summary.\t"+dtf.format(LocalDateTime.now()));
                        executeACAT(output_folder,rscript);
                        get_result(output_folder);
                        System.out.println("Running ACAT with the weights of summary finished.\t"+dtf.format(LocalDateTime.now()));
                    }
                }else {
                    switch (format){
                        case "csv":
                            plink_file=csv2tfile_all(weights_file,genotype_file,phenotype_file,"1,2",output_folder);
                            summary_file=plink2summary(plink,plink_file,"tfile",output_folder);
                            columns="1,3,9";
                            break;
                        case "tfile":
                            summary_file=plink2summary(plink,plink_file,"tfile",output_folder);
                            columns="1,3,9";
                            break;
                        case "bfile":
                            summary_file=plink2summary(plink,plink_file,"bfile",output_folder);
                            columns="1,3,9";
                            break;
                        case "plink":
                            summary_file=plink2summary(plink,plink_file,"file",output_folder);
                            columns="1,3,9";
                            break;
                    }
                    File tmp=new File(output_folder+"acat_tmp/");
                    if(!tmp.exists()){
                        if(!tmp.mkdir()){
                            System.out.println("tmp folder making failed!");
                            System.exit(0);
                        }
                    }

                    BufferedWriter bw=new BufferedWriter(new FileWriter(output_folder+"ACAT_ALL.result"));
                    bw.write("Gene\tPval\n");
                    BufferedReader br=new BufferedReader(new FileReader(genelist));
                    String line;
                    while((line = br.readLine())!= null){
                        try {
                            count=0;
                            Gene_ENSG_ID=line.trim();
                            EN_SNPs_RSID_Dic=new HashMap<>();
                            SNP_EN_Weights_Dic=new HashMap<>();
                            EN_SNPs_Ref_Dic=new HashMap<>();
                            EN_SNPs_Alt_Dic=new HashMap<>();
                            get_en_info(weights_file);
//                            switch (format){
//                                case "csv":
//                                    plink_file=csv2tfile(genotype_file,phenotype_file,"1,2",output_folder);
//                                    summary_file=plink2summary(plink,plink_file,"tfile",output_folder);
//                                    columns="1,3,9";
//                                    break;
//                                case "tfile":
//                                    summary_file=plink2summary(plink,plink_file,"tfile",output_folder);
//                                    columns="1,3,9";
//                                    break;
//                                case "bfile":
//                                    summary_file=plink2summary(plink,plink_file,"bfile",output_folder);
//                                    columns="1,3,9";
//                                    break;
//                                case "plink":
//                                    summary_file=plink2summary(plink,plink_file,"file",output_folder);
//                                    columns="1,3,9";
//                                    break;
//                            }
                            get_weight_p(summary_file,Integer.parseInt(window_size),columns,output_folder+"acat_tmp/");
                            if(count!=0){
                                executeACAT(output_folder+"acat_tmp/",rscript);
                                String line1=get_result1(output_folder+"acat_tmp/");
                                bw.write(line1+"\n");
                            }
                        }catch (Exception e){
                            bw.write(Gene_ENSG_ID+"\tNA\n");
                        }finally {
                            System.out.println(Gene_ENSG_ID);
                        }
                    }
                    br.close();
                    bw.close();
                    try{
//                        File file1 = new File(output_folder+"ACAT.R");
//                        File file2 = new File(output_folder+"ACAT.result");
//                        File file3 = new File(output_folder+"weight_p.txt");
//                        File file4 = new File(output_folder+"plink.assoc");
//                        File file5 = new File(output_folder+"plink.tfam");
//                        File file6 = new File(output_folder+"plink.tped");

//                        if(file1.exists())
//                            file1.delete();
//                        if(file2.exists())
//                            file2.delete();
//                        if(file3.exists())
//                            file3.delete();
//                        if(file4.exists())
//                            file4.delete();
//                        if(file5.exists())
//                            file5.delete();
//                        if(file6.exists())
//                            file6.delete();
                        File[] files = tmp.listFiles();
                        if (files != null) {
                            for (File file : files) {
                                file.delete();
                            }
                        }
                        tmp.delete();

                    }catch(Exception e){
                        e.printStackTrace();
                    }
                    System.out.println("All genes are analyzed successfully!");
                }

                System.out.println("=================END=============");
            }
        }else {
            System.out.println("The function " + function + " doesn't exist.");
        }
    }

}