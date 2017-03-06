/**
 * Project: GNormPlus
 * Function: Species recognition and Species assignment
 */

package GNormPluslib;

import bioc.BioCAnnotation;
import bioc.BioCCollection;
import bioc.BioCDocument;
import bioc.BioCLocation;
import bioc.BioCPassage;

import bioc.io.BioCDocumentWriter;
import bioc.io.BioCFactory;
import bioc.io.woodstox.ConnectorWoodstox;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.BreakIterator;
import java.time.LocalDate;
import java.time.ZoneId;

import javax.xml.stream.XMLStreamException;

import org.tartarus.snowball.SnowballStemmer;
import org.tartarus.snowball.ext.englishStemmer;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

public class SR 
{
	public void SpeciesRecognition(String Filename,String FilenameBioC) throws IOException, XMLStreamException
	{
		/** Recognizing Species Names: SP */
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++) /** PMIDs : i */
		{
			String Pmid = GNormPlus.BioCDocobj.PMIDs.get(i);
			HashMap<String, String> SPID_hash = new HashMap<String, String>();
			PrefixTree PT_Genus = new PrefixTree();
			PrefixTree PT_Strain = new PrefixTree();
			ArrayList<String> tmpStart = new ArrayList<String>();
			ArrayList<String> tmpLast = new ArrayList<String>();
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageName= GNormPlus.BioCDocobj.PassageNames.get(i).get(j); // Passage name
				int PassageOffset = GNormPlus.BioCDocobj.PassageOffsets.get(i).get(j); // Passage offset
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
				
				/** Species recognition */
				ArrayList<String> locations = GNormPlus.PT_Species.SearchMentionLocation(PassageContext,"Species"); /** PT_Species */
				for (int k = 0 ; k < locations.size() ; k++)
				{
					String anno[]=locations.get(k).split("\t");
					int start= Integer.parseInt(anno[0]);
	        		int last= Integer.parseInt(anno[1]);
	        		
	        		// For anti-serum filtering
	        		String ForwardSTR="";
	        		String BackwardSTR="";
	        		if(start>21)
	        		{
	        			ForwardSTR = PassageContext.substring(start-21,last);
	        		}
	        		else
	        		{
	        			ForwardSTR = PassageContext.substring(0,last);
	        		}
	        		if(PassageContext.length()>last+21)
	        		{
	        			BackwardSTR = PassageContext.substring(start,last+21);
	        		}
	        		else
	        		{
	        			BackwardSTR = PassageContext.substring(start,PassageContext.length());
	        		}
	        			        		
	        		String mention = anno[2];
	        		String id = anno[3];
	        		String mention_tmp=mention.toLowerCase();
	        		mention_tmp = mention_tmp.replaceAll("([^A-Za-z0-9@ ])", "\\\\$1");
	        		if(BackwardSTR.toLowerCase().matches(mention_tmp+"[0-9].*")){} // filtered: Bee1p
	        		else if(ForwardSTR.toLowerCase().matches(".*(anti|antibody|antibodies|serum|polyclonal|monoclonal|igg)[\\W\\-\\_]+"+mention_tmp)) {}//filtering : antibody
	        		else if(BackwardSTR.toLowerCase().matches(mention_tmp+"[\\W\\-\\_]+(anti|antibody|antibodies|serum|polyclonal|monoclonal|igg).*")){} //filtering : antibody
	        		else if(BackwardSTR.toLowerCase().matches(mention_tmp+"[\\W\\-\\_]+[A-Za-z0-9]+[\\W\\-\\_]+(anti|antibody|antibodies|serum|polyclonal|monoclonal|igg).*")){} //filtering : antibody
	        		else if(!id.equals("NA"))
	        		{
	        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tSpecies\t"+id); //paragraph
	        			SPID_hash.put(id, "");
	        			tmpStart.add(j+"\t"+start);
	        			tmpLast.add(j+"\t"+last);
	        		}
				}
				
				/** Cell Line recognition */
				locations = GNormPlus.PT_Cell.SearchMentionLocation(PassageContext,"Cell"); /** PT_Cell */
				for (int k = 0 ; k < locations.size() ; k++)
				{
					String anno[]=locations.get(k).split("\t");
					int start= Integer.parseInt(anno[0]);
	        		int last= Integer.parseInt(anno[1]);
	        		String mention = anno[2];
	        		String id = anno[3];
	        		if(!GNormPlus.BioCDocobj.Annotations.get(i).get(j).contains(start+"\t"+last+"\t"+mention+"\tSpecies\t"+id)) //already exists in species result
	        		{
	        			int last40=0;
	        			if(PassageContext.length()>=last+40)
	        			{
	        				last40=last+40;
	        			}
	        			else
	        			{
	        				last40=PassageContext.length();
	        			}
	        			
	        			String patt="[\\W\\-]cell([\\- ]*line|)[s]*[\\W\\-]";
	    				Pattern ptmp = Pattern.compile(patt);
	    				Matcher mtmp = ptmp.matcher(PassageContext.substring(last, last40).toLowerCase());
	    				if(mtmp.find())
	    				{
	    					GNormPlus.BioCDocobj.Annotations.get(i).get(j).add(start+"\t"+last+"\t"+mention+"\tCell\t"+id); //paragraph
	    				}
	        		}
				}
			}
			
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j);
				ArrayList<String> AnnotationInPassage = new ArrayList<String>();
				if(!GNormPlus.BioCDocobj.Annotations.get(i).get(j).isEmpty())
				{
					AnnotationInPassage = GNormPlus.BioCDocobj.Annotations.get(i).get(j);
				}
				
				/** Recognizing Genus (from recognized species names) */
				HashMap<String, String> GenusNames = new HashMap<String, String>();
				for(String ID: SPID_hash.keySet())
				{
					if(GNormPlus.GenusID_hash.containsKey(ID))
					{
						GenusNames.put(ID,GNormPlus.GenusID_hash.get(ID));
					}
				}
				GenusNames.put("3702", "arabidopsis");
				GenusNames.put("4932", "saccharomyces");
				GenusNames.put("562", "escherichia");
				GenusNames.put("7227", "drosophila");
				GenusNames.put("8355", "xenopus");
				
				PT_Genus.Hash2Tree(GenusNames);
				ArrayList<String> locations_Genus = PT_Genus.SearchMentionLocation(PassageContext,"Genus"); /** PT_Genus*/
				for (int k = 0 ; k < locations_Genus.size() ; k++)
				{
					String anno[]=locations_Genus.get(k).split("\t");
					String start= anno[0];
	        		String last= anno[1];
	        		String mention = anno[2];
	        		String id = anno[3];
	        		if(!tmpStart.contains(j+"\t"+start)) //tmpStart
	        		{
	        			AnnotationInPassage.add(start+"\t"+last+"\t"+mention+"\tGenus\t"+id);
	        		}
	        	}
				
				/** Recognizing Strain */
				HashMap<String, String> StrainNames = new HashMap<String, String>();
				for(String ID: SPID_hash.keySet())
				{
					if(GNormPlus.StrainID_hash.containsKey(ID))
					{
						StrainNames.put(ID,GNormPlus.StrainID_hash.get(ID));
					}
				}
				PT_Strain.Hash2Tree(StrainNames);
				ArrayList<String> locations_Strain = PT_Strain.SearchMentionLocation(PassageContext,"Genus"); /** PT_Genus*/
				for (int k = 0 ; k < locations_Strain.size() ; k++)
				{
					String anno[]=locations_Strain.get(k).split("\t");
					String start= anno[0];
	        		String last= anno[1];
	        		String mention = anno[2];
	        		String id = anno[3];
	        		if(!tmpLast.contains(j+"\t"+last)) //tmpLast
	        		{
	        			AnnotationInPassage.add(start+"\t"+last+"\t"+mention+"\tStrain\t"+id);
	        		}
				}
				GNormPlus.BioCDocobj.Annotations.get(i).set(j, AnnotationInPassage);
			}
		}
		//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false); //save in BioC file
	}
	public void SpeciesAssignment(String Filename,String FilenameBioC) throws IOException, XMLStreamException
	{
		BreakIterator iterator = BreakIterator.getSentenceInstance(Locale.US);	
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++) /** PMIDs : i */
		{
			HashMap<String, String> PrefixIDTarget_hash = new HashMap<String, String>();
			PrefixIDTarget_hash.put("9606", "h");
			PrefixIDTarget_hash.put("10090", "m");
			PrefixIDTarget_hash.put("10116", "r");
			PrefixIDTarget_hash.put("4932", "y");
			PrefixIDTarget_hash.put("7227", "d");
			PrefixIDTarget_hash.put("7955", "z|zf|Zf|dr|Dr");
			PrefixIDTarget_hash.put("3702", "at|At");
			
			HashMap<String, Double> SP2Num_hash = new HashMap<String, Double>();
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				ArrayList<String> Annotations_Species = new ArrayList<String>();
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					if(anno.length==5) //Species
	        		{
						if(!PrefixIDTarget_hash.containsKey(anno[4]))
						{
							PrefixIDTarget_hash.put(anno[4],GNormPlus.PrefixID_hash.get(anno[4])); // taxid -> prefix
						}
						if(j == 0)//title
		        		{
		        			if(SP2Num_hash.containsKey(anno[4]))
		        			{
		        				SP2Num_hash.put(anno[4], SP2Num_hash.get(anno[4])+2);
		        			}
		        			else
		        			{
		        				if(GNormPlus.TaxFreq_hash.containsKey(anno[4]))
		        				{
		        					SP2Num_hash.put(anno[4], 2 + GNormPlus.TaxFreq_hash.get(anno[4]));
		        				}
		        				else
		        				{
		        					SP2Num_hash.put(anno[4], 2.0);
		        				}
		        			}
		        			// Virus -> Human
		        			if(GNormPlus.SP_Virus2Human_hash.containsKey(anno[4]))
	        				{
		        				if(SP2Num_hash.containsKey("9606"))
			        			{
			        				SP2Num_hash.put("9606", SP2Num_hash.get("9606")+2);
			        			}
			        			else
			        			{
			        				SP2Num_hash.put("9606", 2 + GNormPlus.TaxFreq_hash.get("9606"));
			        			}
	        				}
		        		}
		        		else
		        		{
		        			if(SP2Num_hash.containsKey(anno[4]))
		        			{
		        				SP2Num_hash.put(anno[4], SP2Num_hash.get(anno[4])+1);
		        			}
		        			else
		        			{
		        				if(GNormPlus.TaxFreq_hash.containsKey(anno[4]))
		        				{
		        					SP2Num_hash.put(anno[4], 1 + GNormPlus.TaxFreq_hash.get(anno[4]));
		        				}
		        				else
		        				{
		        					SP2Num_hash.put(anno[4], 1.0);
		        				}
		        			}
		        			// Virus -> Human
		        			if(GNormPlus.SP_Virus2Human_hash.containsKey(anno[4]))
	        				{
		        				if(SP2Num_hash.containsKey("9606"))
			        			{
			        				SP2Num_hash.put("9606", SP2Num_hash.get("9606")+1);
			        			}
			        			else
			        			{
			        				SP2Num_hash.put("9606", 1 + GNormPlus.TaxFreq_hash.get("9606"));
			        			}
	        				}
		        		}
	        		}
				}
			}
			String MajorSP="9606";
			double MaxSP=0;
			for(String tid : SP2Num_hash.keySet())
			{
				if(SP2Num_hash.get(tid)>MaxSP)
				{
					MajorSP=tid;
					MaxSP=SP2Num_hash.get(tid);
				}
			}
			
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				String PassageContext = GNormPlus.BioCDocobj.PassageContexts.get(i).get(j); // Passage context
				int PassageOffset = GNormPlus.BioCDocobj.PassageOffsets.get(i).get(j); // Passage offset
				iterator.setText(PassageContext);
				ArrayList<Integer> Sentence_offsets = new ArrayList<Integer>();
				int Sent_start = iterator.first();
				for (int Sent_last = iterator.next(); Sent_last != BreakIterator.DONE; Sent_start = Sent_last, Sent_last = iterator.next()) 
				{
					Sentence_offsets.add(Sent_start);
				}
				
				HashMap<Integer,String> Annotations_Gene_hash = new HashMap<Integer,String>();
				ArrayList<String> Annotations_Species = new ArrayList<String>();
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					if(anno.length==5) //Species
	        		{
						Annotations_Species.add(GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k));
	        		}
	        		else //Gene : if(anno.length==3)
	        		{
	        			//String mention = PassageContext.substring(Integer.parseInt(anno[0]), Integer.parseInt(anno[1]));
	        			Annotations_Gene_hash.put(k,GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k)); // k -> Gene Annotation
	        		}
				}

				for (int k : Annotations_Gene_hash.keySet()) // k is the index of GNormPlus.BioCDocobj.Annotations.get(i).get(j) 
    			{
    				String anno[] = Annotations_Gene_hash.get(k).split("\t");
    				int G_Start= Integer.parseInt(anno[0]);
	        		int G_Last= Integer.parseInt(anno[1]);
	        		String G_mentions = anno[2];
	        		String G_type = anno[3];
	        		String G_mention_list[]=G_mentions.split("\\|");
	        		String G_mention=G_mention_list[0]; // only use the first term to detect species ; should be updated after SimConcept
	        		
	        		/** 1. prefix */
	        		boolean SPfound = false;
	        		for(String taxid: PrefixIDTarget_hash.keySet())
	        		{
	        			Pattern ptmp = Pattern.compile("^("+PrefixIDTarget_hash.get(taxid)+")([A-Z].*)$");
						Matcher mtmp = ptmp.matcher(G_mention);
						if(mtmp.find())
						{
							String MentionWoPrefix=mtmp.group(2);
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, anno[0]+"\t"+anno[1]+"\t"+anno[2]+"|"+MentionWoPrefix+"\t"+anno[3]+"\tPrefix:"+taxid);
	        				SPfound=true;
	        				break;
	        			}
	        		}
	        		
	        		/**
	        		 *  2. Co-occurring word
	        		 *  boundary : 
	        		 *  Sentence Start: Sentence_offsets.get(Target_Sentence)
	        		 *  Sentence Last: Sentence_offsets.get(Target_Sentence+1)
	        		 */
	        		//Find the target sentence
	        		int Target_Sentence=0;
	        		if(SPfound == false) // 1. left : Closed to start of the gene mention 
	        		{
	        			for(int s=0;s<Sentence_offsets.size();s++)
	        		
		        		{
		        			int Sentence_last=1000000;
		        			if(s<Sentence_offsets.size()-1)
		        			{
		        				Sentence_last=Sentence_offsets.get(s+1);
		        			}
		        			if(G_Start<Sentence_last)
		        			{
		        				Target_Sentence=s;
		        				break;
		        			}
		        		}
	        		}
	        		int Sentence_Start = Sentence_offsets.get(Target_Sentence);
	        		int Sentence_Last = 1000000;
	        		if(Sentence_offsets.size() > Target_Sentence+1){ Sentence_Last = Sentence_offsets.get(Target_Sentence+1); }
	        		if(SPfound == false) // 1. left : Closed to start of the gene mention 
	        		{
	        			int closet_Sp_Start=0;
	        			for(int sp=0;sp<Annotations_Species.size();sp++) // Find the closet species
		        		{
		        			String AnnoSp[]=Annotations_Species.get(sp).split("\t");
		        			int Sp_Start = Integer.parseInt(AnnoSp[0]);
			        		String taxid = AnnoSp[4];
			        		
			        		if(Sp_Start <= G_Start && Sp_Start >= Sentence_Start && Sp_Start >closet_Sp_Start)
			        		{
			        			closet_Sp_Start=Sp_Start;
			        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tLeft:"+taxid);
		        				SPfound=true;
			        		}
			        	}
		        	}
	        		if(SPfound == false) // 2. right : Closed to last of the gene mention
	        		{
	        			int closet_Sp_Last=1000000;
	        			for(int sp=0;sp<Annotations_Species.size();sp++) // Find the closet species
		        		{
		        			String AnnoSp[]=Annotations_Species.get(sp).split("\t");
		        			int Sp_Last = Integer.parseInt(AnnoSp[1]);
			        		String taxid = AnnoSp[4];
			        		
			        		if(Sp_Last >= G_Last && Sp_Last <= Sentence_Last && Sp_Last < closet_Sp_Last)
			        		{
			        			closet_Sp_Last=Sp_Last;
			        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tRight:"+taxid);
		        				SPfound=true;
			        		}
			        	}
	        		}
	        		
	    			/** 3. Focus species */
	        		if(SPfound == false) // 2. right : Closed to last of the gene mention
	        		{
	        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, Annotations_Gene_hash.get(k)+"\tFocus:"+MajorSP);
        			}
    			}
			}
		}
		//GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false);		
	}
	public void SpeciesAssignment(String Filename,String FilenameBioC,String FocusSpecies) throws IOException, XMLStreamException
	{
		for (int i = 0; i < GNormPlus.BioCDocobj.PMIDs.size(); i++) /** PMIDs : i */
		{
			for (int j = 0; j < GNormPlus.BioCDocobj.PassageNames.get(i).size(); j++) /** Paragraphs : j */
			{
				for (int k = 0; k < GNormPlus.BioCDocobj.Annotations.get(i).get(j).size(); k++) // Annotation : k
				{
					String anno[] = GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k).split("\t");
					if(anno.length==5) //Species
	        		{
						//
	        		}
	        		else //Gene : if(anno.length==3)
	        		{
	        			/** 1. prefix */
		        		boolean SPfound = false;
		        		Pattern ptmp = Pattern.compile("^("+GNormPlus.PrefixID_hash.get(FocusSpecies)+")([A-Z].*)$");
						Matcher mtmp = ptmp.matcher(anno[2]);
						if(mtmp.find())
						{
							String MentionWoPrefix=mtmp.group(2);
							GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k, anno[0]+"\t"+anno[1]+"\t"+anno[2]+"|"+MentionWoPrefix+"\t"+anno[3]+"\tPrefix:"+FocusSpecies);
	        				SPfound=true;
	        			}
		        		if(SPfound == false)
		        		{
		        			GNormPlus.BioCDocobj.Annotations.get(i).get(j).set(k,  GNormPlus.BioCDocobj.Annotations.get(i).get(j).get(k)+"\tFocus:"+FocusSpecies);
		        		}
	        		}
				}
			}
		}
		GNormPlus.BioCDocobj.BioCOutput(Filename,FilenameBioC,GNormPlus.BioCDocobj.Annotations,false);
	}
}