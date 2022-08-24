import os
import sys
import time
import glob
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
from IPython.display import clear_output


class taxaassign :
    def __init__( self, directory ) :
        print( "Pipeline initiating.." )
        global file_vsearch, rdp_classifier
        global dir_user, main_dir, dir_merged, dir_filtered, dir_OTUs, dir_tax_assign, dir_result_summary
        self.directory = directory

        file_vsearch = "vsearch"
        rdp_classifier = "rdp_classifier"

        dir_user = os.path.abspath( os.path.join( os.getcwd() ) )
        main_dir = dir_user + "/TaxaAssignPy/"
        result_dir = main_dir + "result/"
        dir_merged = result_dir + "fastq_merged/"
        dir_filtered = result_dir + "fastq_filtered/"
        dir_OTUs = result_dir + "OTU_define/"
        dir_tax_assign = result_dir + "taxonomy_assign/"
        dir_result_summary = result_dir + "result_summary/"

        if os.path.exists( main_dir ):
            shutil.rmtree( main_dir )
        shutil.copytree( directory, main_dir)

        if os.path.exists( result_dir ):
            shutil.rmtree( result_dir )
        os.makedirs( result_dir )
        os.makedirs( dir_merged )
        os.makedirs( dir_filtered )
        os.makedirs( dir_OTUs )
        os.makedirs( dir_tax_assign )
        os.makedirs( dir_result_summary )
        

    @staticmethod
    def fastq_merge():
        clear_output( wait=True )
        sample_files = main_dir + "*.fastq"
        lst_sample_files = glob.glob( sample_files )
        lst_sample_files.sort()

        if len(lst_sample_files) == 0:
            sample_files = main_dir + "*.fastq.gz"
            lst_sample_files = glob.glob( sample_files )
            lst_sample_files.sort()
            for file in lst_sample_files:
                command = "gzip -d %s" % ( file )
                (exitstatus, outtext) = subprocess.getstatusoutput( command )

        sample_files = main_dir + "*.fastq"
        lst_sample_files = glob.glob( sample_files )
        lst_sample_files.sort()

        for raw_fastq in range( len( lst_sample_files ) ):
            if raw_fastq % 2 == 1:
                continue
            file_f = lst_sample_files[ raw_fastq ]
            file_r = lst_sample_files[ raw_fastq + 1 ]
            file_name = file_f.split( "/" )[ -1 ].split( "_" )[ 0 ]
            merged_out = dir_merged + "%s_merged.fastq" % file_name
            command = "%s --fastq_mergepairs %s --reverse %s --fastqout %s" % ( file_vsearch, file_f, file_r, merged_out )
            print( command )
            ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
            print( outtext )
            print( "\n------------------------------------------------------------------------\n" )

    @staticmethod
    def fastq_filtering():
        clear_output( wait=True )
        merged_samples = os.listdir( dir_merged )
        merged_samples.sort()

        for merged_sample in merged_samples:
            merged_file = dir_merged + merged_sample
            filtered_out = dir_filtered + merged_sample.replace( "_merged", "_filtered" ).replace( ".fastq", ".fasta" )
            command = "%s --fastq_filter %s --fastaout %s" % ( file_vsearch, merged_file, filtered_out )
            print( command )
            ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
            print( outtext )
            print( "\n------------------------------------------------------------------------\n" )


    @staticmethod
    def OTU_defining():
        clear_output( wait=True )
        total_file = dir_OTUs + "total.fasta"
        f = open( total_file, "w")
        filtered_samples = os.listdir( dir_filtered )
        filtered_samples.sort()

        for filtered_sample in filtered_samples:
            filtered_file = dir_filtered + filtered_sample
            file_name = filtered_sample.replace( "_filtered", "" ).replace( ".fasta", "" ).replace( "-", "_" ).replace( " ", "_" )

            for seq in SeqIO.parse( filtered_file, "fasta" ):
                to_write = ">%s[%s]\n%s\n" % ( file_name, seq.id, seq.seq )
                f.write( to_write )

        f.close()

        ###   Dereplication   ###
        derep_file = total_file.replace(".fasta", "_derep.fasta")
        command = "%s --derep_fulllength %s --output %s --sizein --sizeout" % (
                    file_vsearch, total_file, derep_file)
        print( command )
        ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
        print( outtext )
        print("\n------------------------------------------------------------------------\n")

        ###   OTU clutering   ###
        otu_file = total_file.replace(".fasta", "_otus.fasta")
        command = "%s --cluster_unoise %s --centroids %s --id %s --minsize %s --strand both" % (
                    file_vsearch, derep_file, otu_file, 0.97, 2)
        print( command )
        ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
        print( outtext )
        print("\n------------------------------------------------------------------------\n")

        ###   Chimera detection   ###
        nch_file = total_file.replace(".fasta", "_otus_nch.fasta")
        rabel = "sample_OTU_"
        command = "%s --uchime_denovo %s --nonchimeras %s --relabel %s" % (
                    file_vsearch, otu_file, nch_file, rabel)
        print( command )
        ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
        print( outtext )
        print( "\n------------------------------------------------------------------------\n" )

        ###   OTU table making   ###
        otu_table = total_file.replace(".fasta", "_otu_table.txt")
        map_uc = total_file.replace(".fasta", "_map.uc")
        command = "%s --usearch_global %s --db %s --otutabout %s --uc %s --id %s --strand both" % (
                          file_vsearch, total_file, nch_file, otu_table, map_uc, 0.97)
        print( command )
        ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
        print( outtext )
        print( "\n------------------------------------------------------------------------\n" )

        nch_file_cp = nch_file.replace( "OTU_define", "result_summary" )
        if os.path.exists( nch_file_cp ):
            shutil.rmtree( nch_file_cp )
        shutil.copy( nch_file, nch_file_cp )
        
        out_table_cp = otu_table.replace( "OTU_define", "result_summary" )
        if os.path.exists( out_table_cp ):
            shutil.rmtree( out_table_cp )
        shutil.copy( otu_table, out_table_cp )

    @staticmethod
    def taxonomy_assigning():
        clear_output( wait=True )
        classifier_in = dir_OTUs + "total_otus_nch.fasta"
        classifier_out = dir_tax_assign + "RDP-classifier_result.txt"
        command = "%s -q %s -o %s" % ( rdp_classifier, classifier_in, classifier_out )
        print( command )
        ( exitstatus, outtext ) = subprocess.getstatusoutput( command )
        print( outtext )

    @staticmethod
    def analysis_statistic():
        clear_output( wait=True )
        classifier_out = dir_tax_assign + "RDP-classifier_result.txt"
        otu_table = dir_OTUs + "total_otu_table.txt"
        filtered_samples = os.listdir( dir_filtered )
        filtered_samples.sort()

        tax_num = 6
        tax_lv_lst = [ "domain", "phylum", "class", "order", "family", "genus" ]
        dic_tax_index = { "domain":0, "phylum":1, "class":2, "order":3, "family":4, "genus":5 }

        ##  Calculating taxonomy results  ##
        print( "Doing statistics...... ", end="" )
        dic_otu_tax = {}
        for num in range( tax_num ):
            dic_otu_tax[ num ] = {}

        f = open( classifier_out, "r" )
        lines = f.readlines()
        f.close()

        f = open( classifier_out.replace( "taxonomy_assign", "result_summary" ), "w")
        for line in lines:
            if line == "":
                continue
            temp_spl = line.replace( "\r", "" ).replace( "\n", "" ).replace( "\"", "" )
            spl_t = temp_spl.split( "\t" )

            for num in range( 4 ):
                if "" in spl_t:
                    spl_t.remove( "" )
                elif "-" in spl_t:
                    spl_t.remove( "-" )

            temp_id = spl_t[ 0 ].split( ";" )[ 0 ]

            f.write( "%s\t" % spl_t[ 0 ] )

            for tax in tax_lv_lst:
                tax_lv = dic_tax_index[ tax ]
                tax_name = "unclassified"
                tax_confi = 0
                if tax in spl_t:
                    tax_index = spl_t.index( tax )
                    tax_name = spl_t[ tax_index - 1 ]
                    tax_confi = spl_t[ tax_index + 1 ]
                if len( tax_name.split( "." ) ) == 3:
                    tax_name = "unclassified"
                elif tax_name.isdigit():
                    tax_name == "unclassified"
                tax_name_w = tax_name
                f.write( "D_%s__%s(%s)" % ( tax_lv, tax_name_w, tax_confi ) )
                if tax_lv != tax_num - 1:
                    f.write( ";" )
                if float(tax_confi) < 0.8:
                    tax_name = "under_confi_lv"
                tax_label = "D_%s__%s" % ( tax_lv, tax_name )
                if tax_label not in dic_otu_tax[ tax_lv ].keys():
                    dic_otu_tax[ tax_lv ][ tax_label ] = [ temp_id ]
                elif tax_label in dic_otu_tax[ tax_lv ].keys():
                    dic_otu_tax[ tax_lv ][ tax_label ].append( temp_id )
            f.write( "\r\n" )
        f.close()
        print( " Done.")

        f = open( otu_table, "r" )
        r = f.read()
        f.close()
        lines = r.split( "\n" )

        dic_otu_order = {}
        dic_sample_order = {}
        matrix = []
        temp_num = 0
        for line in lines:
            if line == "":
                continue
            spl_ = line.split( "\t" )
            if line.startswith( "#OTU" ):
                for name in spl_:
                    dic_sample_order[ name ] = spl_.index( name )
                continue
            matrix.append( spl_[ 1: ] )
            dic_otu_order[ spl_[ 0 ] ] = temp_num
            temp_num += 1

        ##  Writing assignment statistic with descending order depend on taxonomy lv  ##
        sample_order = list( dic_sample_order.keys() )
        sample_order.sort()
        sample_order = sample_order[1:]
        
        for num in dic_otu_tax.keys():
            print( "%s taxonomy level(%s) file writing......" % ( num, tax_lv_lst[ num ] ), end="" )
            summary_matrix = []
            summary_temp_matrix = [ "File_ID" ]
            for k in dic_otu_tax[ num ].keys():
                summary_temp_matrix.append( k )
            summary_matrix.append( summary_temp_matrix )

            dic_temp_num = {}
            count_matrix = []
            for f_name in sample_order:
                summary_temp_matrix = [ f_name ]

                dic_temp_num[f_name] = 0
                temp_matrix = []
                for k, v in dic_otu_tax[ num ].items():
                    temp_num = 0
                    for z in v:
                        if z not in dic_otu_order.keys():
                            continue
                        temp_num += int( matrix[ dic_otu_order[ z ] ][ dic_sample_order[ f_name ] - 1 ] )
                    summary_temp_matrix.append( temp_num )
                    temp_matrix.append( temp_num )
                    dic_temp_num[ f_name ] += temp_num
                summary_matrix.append( summary_temp_matrix )
                count_matrix.append( temp_matrix )

            summary_temp_matrix = [ "File_ID" ]
            for k in dic_otu_tax[ num ].keys():
                summary_temp_matrix.append( k )
            summary_matrix.append( summary_temp_matrix )

            for f_name in sample_order:
                summary_temp_matrix = [ f_name ]
                for i in count_matrix:
                    if count_matrix.index( i ) != int( dic_sample_order[ f_name ] ) - 1:
                        continue
                    for j in i:
                        count_num = float( j )
                        summary_temp_matrix.append( str( ( count_num / dic_temp_num[ f_name ] ) ) )
                    summary_matrix.append( summary_temp_matrix )

            df = pd.DataFrame( summary_matrix )

            total_sample_num = len( filtered_samples )

            total_sum_matrix = [ 0 for i in range( len( df.columns ) ) ]
            for i in range( 1, total_sample_num + 1 ):
                for j in range( 1, len( df.columns ) ):
                    total_sum_matrix[ j ] += int( df.iloc[ i, j ] )

            total_sum_matrix[ 0 ] = sum( total_sum_matrix[ 1: ] )
            df.loc[ total_sample_num + 1 ] = total_sum_matrix
            df = df.sort_values( by = total_sample_num + 1, axis=1, ascending=False )

            total_sum_str = []
            for x in list( df.loc[ total_sample_num + 1 ] ):
                total_sum_str.append( str( x ) )
            re_from = "\t".join( total_sum_str )
            re_to = "\t".join( list( df.loc[ 0 ] ) )

            to_write = df.to_csv( sep="\t", index = False ).replace( re_from, "%s\n\n\n\n%s" % (
                    re_from.replace( str( total_sum_matrix[ 0 ] ) + "\t", "Total\t" ), re_to))
            to_write = "\r\n".join( to_write.split( "\n" )[ 1: ] )
            f = open( dir_result_summary + "%s_rdp-classifier.txt" % num, "w" )
            f.write( to_write )
            f.close()
            print( " Done." )

    @staticmethod
    def result_summary():
        clear_output( wait=True )
        total_file = dir_OTUs + "total.fasta"
        nch_file = dir_OTUs + "total_otus_nch.fasta"
        filtered_samples = os.listdir( dir_filtered )
        filtered_samples.sort()
        sample_files = main_dir + "*.fastq"
        raw_files = glob.glob( sample_files )
        raw_files.sort()

        tax_num = 6
        tax_lv_lst = [ "domain", "phylum", "class", "order", "family", "genus" ]

        ##  Summaey with text  ##
        total_sample_num = len( filtered_samples )

        total_read_num = 0
        for file in raw_files:
            for seq in SeqIO.parse( file, "fastq" ):
                total_read_num += 1
        total_read_num = int( total_read_num / 2 )

        filtered_read_num = 0
        for seq in SeqIO.parse( total_file, "fasta" ):
            filtered_read_num += 1

        total_OTU_num = 0
        for seq in SeqIO.parse( nch_file, "fasta" ):
            total_OTU_num += 1

        total_taxo_num = 0
        for num in range( tax_num ):
            f = open( dir_result_summary + "%s_rdp-classifier.txt" % num, "r" )
            lines = f.readlines()
            f.close()
            for line in lines:
                spl_ = line.split( "\t" )
                temp_taxo_num = len( spl_ ) - 1
                if "D_%s__under_confi_lv" % num in spl_:
                    temp_taxo_num -= 1
                break
            if total_taxo_num == 0:
                total_taxo_num = temp_taxo_num
            elif total_taxo_num < temp_taxo_num:
                total_taxo_num = temp_taxo_num

        f = open( dir_result_summary + "overall_results_summary.txt", "w" )
        to_write = "Result Summary\n%-20s:\t%s\n%-20s:\t%s\n%-20s:\t%s\n%-20s:\t%s\n%-20s:\t%s\n" % (
                    "Total sample #", total_sample_num, "Total read #", total_read_num, "Filtered read #", filtered_read_num, 
                    "Defined OTU #", total_OTU_num, "Assigned taxonomy #", total_taxo_num)
        f.write( to_write )
        f.close()

        print ( to_write )
        print ( "---------------------------------------------------------------------" )

        ##  Summaey with figure  ##
        for num in range( tax_num ):
            f = open( dir_result_summary + "%s_rdp-classifier.txt" % num, "r" )
            r = f.read()
            f.close()
            count_table = r.split( "\n" * 4 )[ 0 ]
            lines = count_table.split( "\n" )
            tax_name_lst = []
            tax_num_lst = []
            for line in lines:
                spl_ = line.split( "\t" )
                if "File_ID" in spl_:
                    tax_name_lst = spl_[ 1: ]
                    continue
                if tax_num_lst == []:
                    for z in spl_[ 1: ]:
                        tax_num_lst.append( int( z ) )
                else:
                    temp_0 = 0
                    for z in spl_[ 1: ]:
                        tax_num_lst[ temp_0 ] += int( z )
                        temp_0 += 1

            print( "\n\n%s level taxonomy" % tax_lv_lst[ num ] )
            data = pd.DataFrame( { "tax_name" : tax_name_lst, "tax_num": tax_num_lst } )
            data = data.sort_values( [ "tax_num" ], ascending=[ True ] )
            data = data.reset_index( drop=True )

            y_values = list( data[ "tax_name" ] )
            x_values = list( data[ "tax_num" ] )
            #plt.style.use( [ 'dark_background' ] )
            fig = plt.figure( figsize=( 15, len( y_values ) / 4 ) )
            ax = fig.add_subplot( 111 )
            ypos = range( len( x_values ) )
            rects = plt.barh( ypos, x_values, align='center', height=0.5 )
            for i, rect in enumerate( rects ):
                ax.text( 1.01 * rect.get_width(), rect.get_y() + rect.get_height() / 2.0, str( x_values[ i ] ), ha='left', va='center')
            plt.yticks( ypos, y_values )
            plt.title( '%s level taxonomy' % tax_lv_lst[ num ] )
            plt.xlabel( 'Number' )
            plt.ylabel( 'Taxonomy' )
            plt.savefig( dir_result_summary + "%s_taxonomy_assign.svg" % num )
            plt.show()

class autorun_all :
    def __init__( self, directory ) :
        print( "Pipeline initiating.." )
        global file_vsearch, rdp_classifier, start_time
        global dir_user, main_dir, dir_merged, dir_filtered, dir_OTUs, dir_tax_assign, dir_result_summary
        self.directory = directory

        file_vsearch = "vsearch"
        rdp_classifier = "rdp_classifier"
        start_time = time.time()

        dir_user = os.path.abspath( os.path.join( os.getcwd() ) )
        main_dir = dir_user + "/TaxaAssignPy/"
        result_dir = main_dir + "result/"
        dir_merged = result_dir + "fastq_merged/"
        dir_filtered = result_dir + "fastq_filtered/"
        dir_OTUs = result_dir + "OTU_define/"
        dir_tax_assign = result_dir + "taxonomy_assign/"
        dir_result_summary = result_dir + "result_summary/"

        if os.path.exists( main_dir ):
            shutil.rmtree( main_dir )
        shutil.copytree( directory, main_dir)

        if os.path.exists( result_dir ):
            shutil.rmtree( result_dir )
        os.makedirs( result_dir )
        os.makedirs( dir_merged )
        os.makedirs( dir_filtered )
        os.makedirs( dir_OTUs )
        os.makedirs( dir_tax_assign )
        os.makedirs( dir_result_summary )

        taxaassign.fastq_merge()
        taxaassign.fastq_filtering()
        taxaassign.OTU_defining()
        taxaassign.taxonomy_assigning()
        taxaassign.analysis_statistic()
        taxaassign.result_summary()

        end_time = time.time() - start_time
        print( "Pipeline time consume : %s Min %s sec" % ( end_time // 60, int( end_time ) % 60 ) )
