/**
 * Project: PubTaggers
 * Function: Main
 */
//
package GNormPluslib;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import GNormPluslib.PrefixTree;
import GNormPluslib.GNR;
import GNormPluslib.SR;

public class GNormPlus
{
	public static BioCDoc BioCDocobj = new BioCDoc();
	public static PrefixTree PT_Species = new PrefixTree();
	public static PrefixTree PT_Cell = new PrefixTree();
	public static PrefixTree PT_CTDGene = new PrefixTree();
	public static PrefixTree PT_Gene = new PrefixTree();
	public static PrefixTree PT_GeneChromosome = new PrefixTree();
	public static HashMap<String, String> ent_hash = new HashMap<String, String>();
	public static HashMap<String, String> GenusID_hash = new HashMap<String, String>();
	public static HashMap<String, String> StrainID_hash = new HashMap<String, String>();
	public static HashMap<String, String> PrefixID_hash = new HashMap<String, String>();
	public static HashMap<String, Double> TaxFreq_hash = new HashMap<String, Double>();
	public static HashMap<String, String> GeneScoring_hash = new HashMap<String, String>();
	public static HashMap<String, Double> GeneScoringDF_hash = new HashMap<String, Double>();
	public static HashMap<String, String> GeneIDs_hash = new HashMap<String, String>();
	public static HashMap<String, String> Normalization2Protein_hash = new HashMap<String, String>();
	public static ArrayList<String> SuffixTranslationMap = new ArrayList<String>();
	public static HashMap<String, String> Pmid2Abb_hash = new HashMap<String, String>();
	public static HashMap<String, String> PmidAbb2LF_lc_hash = new HashMap<String, String>();
	public static HashMap<String, String> PmidLF2Abb_lc_hash = new HashMap<String, String>();
	public static HashMap<String, String> PmidAbb2LF_hash = new HashMap<String, String>();
	public static HashMap<String, String> PmidLF2Abb_hash = new HashMap<String, String>();
	public static HashMap<String, String> Pmid2ChromosomeGene_hash = new HashMap<String, String>();
	public static HashMap<String, String> SimConceptMention2Type_hash = new HashMap<String, String>();
	public static HashMap<String, String> Filtering_hash = new HashMap<String, String>();
	public static HashMap<String, String> SP_Virus2Human_hash = new HashMap<String, String>();
	public static void main(String [] args) throws IOException, InterruptedException, XMLStreamException 
	{
		String InputFolder="input";
		String OutputFolder="output";
		String SetupFile = "setup.txt";
		String FocusSpecies = "";
		if(args.length<2)
		{
			System.out.println("\n$ java -Xmx10G -Xms10G -jar GNormPlus.jar [InputFolder] [OutputFolder] [SetupFile]");
			System.out.println("[InputFolder] Default : input");
			System.out.println("[OutputFolder] Default : output");
			System.out.println("[SetupFile] Default : setup.txt\n\n");
		}
		else
		{
			/*
			 * Parameters
			 */
			InputFolder=args[0];
			OutputFolder=args[1];
			if(args.length>2)
			{
				SetupFile = args[2];
			}
			if(args.length>=4)
			{
				FocusSpecies=args[3];
			}
		}
		
		HashMap<String, String> setup_hash = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(SetupFile));
		String line="";
		Pattern ptmp = Pattern.compile("^	([A-Za-z0-9]+) = ([^ \\t\\n\\r]+)$");
		while ((line = br.readLine()) != null)  
		{
			Matcher mtmp = ptmp.matcher(line);
			if(mtmp.find())
			{
				setup_hash.put(mtmp.group(1), mtmp.group(2));
			}
		}
		br.close();
		if(!setup_hash.containsKey("GeneIDMatch"))
		{
			setup_hash.put("GeneIDMatch","True");
		}
		
		/*
		 * Time stamp - start : All
		 */
		double startTime,endTime,totTime;
		startTime = System.currentTimeMillis();//start time
	
		/* 
		 * Start & Load Dictionary
		 */
		String TrainTest = "Test";
		if(setup_hash.containsKey("TrainTest"))
		{
			TrainTest = setup_hash.get("TrainTest");
		}
		
		System.out.print("Loading Gene Dictionary : Processing ... \r");
		/** Load Dictionary */
		{
			/** CTDGene  */			
			PT_CTDGene.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_CTDGene.txt");
			
			/** ent */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/ent.rev.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t"); //&#x00391;	Alpha
				ent_hash.put(l[0], l[1]);
			}
			br.close();	

			/** Species */
			PT_Species.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_Species.txt");
			
			/** Cell */
			PT_Cell.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_Cell.txt");
			
			/** Genus */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/SPGenus.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				GenusID_hash.put(l[0], l[1]); // tax id -> Genus
			}
			br.close();	
			
			/** Strain */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/SPStrain.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				StrainID_hash.put(l[0], l[1]); // tax id -> strain
			}
			br.close();
			/** Prefix */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/SPPrefix.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				PrefixID_hash.put(l[0], l[1]); //tax id -> prefix
			}
			br.close();
			PrefixID_hash.put("9606", "h");
			PrefixID_hash.put("10090", "m");
			PrefixID_hash.put("10116", "r");
			PrefixID_hash.put("4932", "y");
			PrefixID_hash.put("7227", "d");
			PrefixID_hash.put("7955", "z|dr|Dr|Zf|zf");
			PrefixID_hash.put("3702", "at|At");
			
			/** Frequency */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/taxonomy_freq.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				TaxFreq_hash.put(l[0], Double.parseDouble(l[1])/200000000); //tax id -> prefix
			}
			br.close();	
			
			/** SP_Virus2Human_hash */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/SP_Virus2HumanList.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				SP_Virus2Human_hash.put(line,"9606");
			}
			br.close();	

			/** SimConcept.MentionType */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/SimConcept.MentionType.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				SimConceptMention2Type_hash.put(l[0], l[1]);
			}
			br.close();	
			
			/** Filtering */
			br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/Filtering.txt"));
			line="";
			while ((line = br.readLine()) != null)  
			{
				Filtering_hash.put(line, "");
			}
			br.close();	
				
			/** Gene */
			if(setup_hash.containsKey("FocusSpecies") && !setup_hash.get("FocusSpecies").equals("All"))
			{
				PT_Gene.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_Gene."+setup_hash.get("FocusSpecies")+".txt");
			}
			else if((!FocusSpecies.equals("")) && (!FocusSpecies.equals("All")))
			{
				PT_Gene.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_Gene."+setup_hash.get("FocusSpecies")+".txt");
			}
			else
			{
				PT_Gene.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_Gene.txt");
			}
			
			
			/** GeneScoring */
			String FileName=setup_hash.get("DictionaryFolder")+"/GeneScoring.txt";
			
			if(setup_hash.containsKey("FocusSpecies") && !setup_hash.get("FocusSpecies").equals("All"))
			{
				FileName = setup_hash.get("DictionaryFolder")+"/GeneScoring."+setup_hash.get("FocusSpecies")+".txt";
			}
			else if((!FocusSpecies.equals("")) && (!FocusSpecies.equals("All")))
			{
				FileName = setup_hash.get("DictionaryFolder")+"/GeneScoring."+FocusSpecies+".txt";
			}
			br = new BufferedReader(new FileReader(FileName));
			line="";
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				GeneScoring_hash.put(l[0], l[1]+"\t"+l[2]+"\t"+l[3]+"\t"+l[4]+"\t"+l[5]+"\t"+l[6]);
			}
			br.close();	
			
			/** GeneScoring.DF */
			FileName=setup_hash.get("DictionaryFolder")+"/GeneScoring.DF.txt";
			if(setup_hash.containsKey("FocusSpecies") && !setup_hash.get("FocusSpecies").equals("All"))
			{
				FileName = setup_hash.get("DictionaryFolder")+"/GeneScoring.DF."+setup_hash.get("FocusSpecies")+".txt";
			}
			else if((!FocusSpecies.equals("")) && (!FocusSpecies.equals("All")))
			{
				FileName = setup_hash.get("DictionaryFolder")+"/GeneScoring.DF."+FocusSpecies+".txt";
			}
			br = new BufferedReader(new FileReader(FileName));
			double Sum = Double.parseDouble(br.readLine());
			while ((line = br.readLine()) != null)  
			{
				String l[]=line.split("\t");
				// token -> idf
				GeneScoringDF_hash.put(l[0], Math.log10(Sum/Double.parseDouble(l[1])));
			}
			br.close();
			
			/** Suffix Translation */
			SuffixTranslationMap.add("alpha-a");
			SuffixTranslationMap.add("alpha-1");
			SuffixTranslationMap.add("a-alpha");
			SuffixTranslationMap.add("a-1");
			SuffixTranslationMap.add("1-alpha");
			SuffixTranslationMap.add("1-a");
			SuffixTranslationMap.add("beta-b");
			SuffixTranslationMap.add("beta-2");
			SuffixTranslationMap.add("b-beta");
			SuffixTranslationMap.add("b-2");
			SuffixTranslationMap.add("2-beta");
			SuffixTranslationMap.add("2-b");
			SuffixTranslationMap.add("gamma-g");
			SuffixTranslationMap.add("gamma-y");
			SuffixTranslationMap.add("g-gamma");
			SuffixTranslationMap.add("y-gamma");
			
			/** GeneID */
			if(setup_hash.containsKey("GeneIDMatch") && setup_hash.get("GeneIDMatch").toLowerCase().equals("true"))
			{
				br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/GeneIDs.txt"));
				line="";
				while ((line = br.readLine()) != null)  
				{
					String l[]=line.split("\t");
					GeneIDs_hash.put(l[0],l[1]);
				}
				br.close();
			}
			
			/** Normalization2Protein */
			if(setup_hash.containsKey("Normalization2Protein") && setup_hash.get("Normalization2Protein").toLowerCase().equals("true"))
			{
				br = new BufferedReader(new FileReader(setup_hash.get("DictionaryFolder")+"/Gene2Protein.txt"));
				line="";
				while ((line = br.readLine()) != null)  
				{
					String l[]=line.split("\t");
					Normalization2Protein_hash.put(l[0],l[1]);
				}
				br.close();
			}
			
			/** GeneChromosome */
			//PT_GeneChromosome.TreeFile2Tree(setup_hash.get("DictionaryFolder")+"/PT_GeneChromosome.txt");
		}
		endTime = System.currentTimeMillis();
		totTime = endTime - startTime;
		System.out.println("Loading Gene Dictionary : Processing Time:"+totTime/1000+"sec");
		
		File folder = new File(InputFolder);
		File[] listOfFiles = folder.listFiles();
		for (int i = 0; i < listOfFiles.length; i++)
		{
			if (listOfFiles[i].isFile()) 
			{
				String InputFile = listOfFiles[i].getName();
				File f = new File(OutputFolder+"/"+InputFile);
				if(f.exists() && !f.isDirectory()) 
				{ 
					System.out.println(InputFolder+"/"+InputFile+" - Done. (The output file exists in output folder)");
				}
				else
				{
					BioCDocobj = new BioCDoc();
					
					/*
					 * Format Check 
					 */
					String Format = "";
					String checkR = BioCDocobj.BioCFormatCheck(InputFolder+"/"+InputFile);
					if(checkR.equals("BioC"))
					{
						Format = "BioC";
					}
					else if(checkR.equals("PubTator"))
					{
						Format = "PubTator";
					}
					else
					{
						System.out.println(checkR);
						System.exit(0);
					}
					
					System.out.print(InputFolder+"/"+InputFile+" - ("+Format+" format) : Processing ... \r");
					/*
					 * GNR
					 */
					
					{
						GNR GNRobj = new GNR();
						
						if(setup_hash.containsKey("IgnoreNER") && setup_hash.get("IgnoreNER").equals("True"))
						{
							if(Format.equals("PubTator"))
							{
								BioCDocobj.PubTator2BioC(InputFolder+"/"+InputFile,"tmp/"+InputFile);
								GNRobj.LoadInputFile("tmp/"+InputFile,"tmp/"+InputFile+".Abb","Train");
							}
							else if(Format.equals("BioC"))
							{
								GNRobj.LoadInputFile(InputFolder+"/"+InputFile,"tmp/"+InputFile+".Abb","Train");
							}
						}
						else
						{
							if(Format.equals("PubTator"))
							{
								BioCDocobj.PubTator2BioC(InputFolder+"/"+InputFile,"tmp/"+InputFile);
								GNRobj.LoadInputFile("tmp/"+InputFile,"tmp/"+InputFile+".Abb",TrainTest);
							}
							else if(Format.equals("BioC"))
							{
								GNRobj.LoadInputFile(InputFolder+"/"+InputFile,"tmp/"+InputFile+".Abb",TrainTest);
							}
							
							GNRobj.FeatureExtraction("tmp/"+InputFile+".data","tmp/"+InputFile+".loca",TrainTest);
							
							if(TrainTest.equals("Test"))
							{
								GNRobj.CRF_test(setup_hash.get("GNRModel"),"tmp/"+InputFile+".data","tmp/"+InputFile+".output","top3"); //top3
								
								if(Format.equals("PubTator"))
								{
									GNRobj.ReadCRFresult("tmp/"+InputFile,"tmp/"+InputFile+".loca","tmp/"+InputFile+".output","tmp/"+InputFile+".GNR.xml",0.005,0.05); //0.005,0.05
									GNRobj.PostProcessing("tmp/"+InputFile,"tmp/"+InputFile+".GNR.xml");
								}
								else if(Format.equals("BioC"))
								{
									GNRobj.ReadCRFresult(InputFolder+"/"+InputFile,"tmp/"+InputFile+".loca","tmp/"+InputFile+".output","tmp/"+InputFile+".GNR.xml",0.005,0.05); //0.005,0.05
									GNRobj.PostProcessing(InputFolder+"/"+InputFile,"tmp/"+InputFile+".GNR.xml");
								}
							}
							else if(TrainTest.equals("Train"))
							{
								//GNRobj.CRF_learn(setup_hash.get("GNRModel"),"tmp/"+InputFile+".data");
							}
						}
					}
					
					/*
					 * SR & SA
					 */
					if(TrainTest.equals("Test"))
					{
						SR SRobj = new SR();
						SRobj.SpeciesRecognition(InputFolder+"/"+InputFile,"tmp/"+InputFile+".SA.xml");
						
						if((!FocusSpecies.equals("")) && (!FocusSpecies.equals("All")))
						{
							if(Format.equals("PubTator"))
							{
								SRobj.SpeciesAssignment("tmp/"+InputFile,"tmp/"+InputFile+".SA.xml",FocusSpecies);
							}
							else if(Format.equals("BioC"))
							{
								SRobj.SpeciesAssignment(InputFolder+"/"+InputFile,"tmp/"+InputFile+".SA.xml",FocusSpecies);
							}
						}
						else if(setup_hash.containsKey("FocusSpecies") && !setup_hash.get("FocusSpecies").equals("All"))
						{
							if(Format.equals("PubTator"))
							{
								SRobj.SpeciesAssignment("tmp/"+InputFile,"tmp/"+InputFile+".SA.xml",setup_hash.get("FocusSpecies"));
							}
							else if(Format.equals("BioC"))
							{
								SRobj.SpeciesAssignment(InputFolder+"/"+InputFile,"tmp/"+InputFile+".SA.xml",setup_hash.get("FocusSpecies"));
							}
						}
						else
						{
							if(Format.equals("PubTator"))
							{
								SRobj.SpeciesAssignment("tmp/"+InputFile,"tmp/"+InputFile+".SA.xml");
							}
							else if(Format.equals("BioC"))
							{
								SRobj.SpeciesAssignment(InputFolder+"/"+InputFile,"tmp/"+InputFile+".SA.xml");
							}
						}
					}
					
					/*
					 * SimConcept
					 */
					{
						SimConcept SCobj = new SimConcept();
						if(TrainTest.equals("TrainSC"))
						{
							SCobj.FeatureExtraction_Train("tmp/"+InputFile+".SC.data");
							SCobj.CRF_learn(setup_hash.get("SCModel"),"tmp/"+InputFile+".SC.data");
						}
						else if(TrainTest.equals("Test"))
						{
							SCobj.FeatureExtraction_Test("tmp/"+InputFile+".SC.data");
							SCobj.CRF_test(setup_hash.get("SCModel"),"tmp/"+InputFile+".SC.data","tmp/"+InputFile+".SC.output");
							
							if(Format.equals("PubTator"))
							{
								SCobj.ReadCRFresult("tmp/"+InputFile,"tmp/"+InputFile+".SC.output","tmp/"+InputFile+".SC.xml");
							}
							else
							{
								SCobj.ReadCRFresult(InputFolder+"/"+InputFile,"tmp/"+InputFile+".SC.output","tmp/"+InputFile+".SC.xml");
							}
							
						}
					}
					
					/*
					 * GN
					 */
					if(TrainTest.equals("Test"))
					{
						GN GNobj = new GN();
						GNobj.PreProcessing4GN(InputFolder+"/"+InputFile,"tmp/"+InputFile+".PreProcessing4GN.xml");
						GNobj.ChromosomeRecognition(InputFolder+"/"+InputFile,"tmp/"+InputFile+".GN.xml");
						if(setup_hash.containsKey("GeneIDMatch") && setup_hash.get("GeneIDMatch").equals("True"))
						{
							
							if(Format.equals("PubTator"))
							{
								GNobj.GeneNormalization("tmp/"+InputFile,"tmp/"+InputFile+".GN.xml",true);
								GNobj.GeneIDRecognition("tmp/"+InputFile,"tmp/"+InputFile+".GN.xml");
								BioCDocobj.BioC2PubTator("tmp/"+InputFile+".GN.xml",OutputFolder+"/"+InputFile);
							}
							else if(Format.equals("BioC"))
							{
								GNobj.GeneNormalization(InputFolder+"/"+InputFile,"tmp/"+InputFile+".GN.xml",true);
								GNobj.GeneIDRecognition(InputFolder+"/"+InputFile,OutputFolder+"/"+InputFile);
							}
						}
						else
						{
							if(Format.equals("PubTator"))
							{
								GNobj.GeneNormalization("tmp/"+InputFile,"tmp/"+InputFile+".GN.xml",false);
								BioCDocobj.BioC2PubTator("tmp/"+InputFile+".GN.xml",OutputFolder+"/"+InputFile);
							}
							else if(Format.equals("BioC"))
							{
								GNobj.GeneNormalization(InputFolder+"/"+InputFile,OutputFolder+"/"+InputFile,false);
							}
						}
					}
					
					/*
					 * remove tmp files
					 */
					if((!setup_hash.containsKey("DeleteTmp")) || setup_hash.get("DeleteTmp").toLowerCase().equals("true"))
					{
						String path="tmp"; 
				        File file = new File(path);
				        File[] files = file.listFiles(); 
				        for (File ftmp:files) 
				        {
				        	if (ftmp.isFile() && ftmp.exists()) 
				            {
				        		if(ftmp.toString().matches("tmp."+InputFile+".*"))
					        	{
				        			ftmp.delete();
					        	}
				        	}
				        }
					}
					
					/*
					 * Time stamp - last
					 */
					endTime = System.currentTimeMillis();
					totTime = endTime - startTime;
					System.out.println(InputFolder+"/"+InputFile+" - ("+Format+" format) : Processing Time:"+totTime/1000+"sec");
				}
			}
		}
	}
}
