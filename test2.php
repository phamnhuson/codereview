<?php
	namespace App\Http\Controllers;
	
	use Illuminate\Http\Request;
	use Validator;
	use Parser;
	use DB;
	
	class BlastController extends Controller
	{
	
		function index()
		{
			return view('blast');
		}
	
		function blast(Request $request)
		{
			$validator = Validator::make($request->all(), [
				'sequence'	=>	'Regex:/^([ATGCatgc\-\n]+)$/'
			]);
			
			if ($validator->fails()) {
			
				return \Redirect::back()
						->withErrors($validator)
						->withInput();
						
			} else {
				$viewData = array();
			
				$viewData['blastResult'] = $this->blastProccess($request->all());
				$viewData['oldInput'] = $request->all();
				return view('blast', $viewData);
				
			}
		}
	
		function blastProccess($request)
		{
			
			$blastTool = storage_path()."/blast+/bin/".$request['tool'];
			$threadsHold = $request['threadshold'];
			$wordSize = $request['wordsize'];
			$targetSeqs = $request['targetseqs'];
			$scores	= $request['scores'];
			
			if (isset($request['sequence'])) {
				$sequence = $request['sequence'];
				$inputFile = storage_path().'/linux/blast.input';
				file_put_contents($inputFile, $sequence);
			}	
			
			$blastDb = '/linux/fasta/nucleotide_db';
			if (isset($request['allowmore'])) {
				$tmpName = md5(uniqid());
				$blastDb = '/linux/fasta-blast/'.$tmpName;
				file_put_contents(storage_path().$blastDb, $request['moresequence']);
				
				// Call make database function
				chdir(storage_path().'/linux/fasta-blast/');
				exec(storage_path().config('app.blast_tool_path')." -in ".$tmpName." -input_type fasta -dbtype nucl -out ".$tmpName);
			}
			
			
			$command = "$blastTool -query $inputFile -db ".storage_path().$blastDb." -outfmt 5";
			
			if ($threadsHold) {
				$command .= " -num-threads $threadsHold";
			}
			
			if ($wordSize) {
				$command .= " -word_size $wordSize";
			}
			
			if ($targetSeqs) {
				$command .= " -max_target_seqs $targetSeqs";
			}
			
			ob_start();
			system($command, $xml);
			$xml = ob_get_clean();
			
			$result = Parser::xml($xml);
			
			$blastResult = array();
			
			if (isset($result['BlastOutput_iterations']) && !empty($result['BlastOutput_iterations']['Iteration']['Iteration_hits'])) {
				
				if (isset($result['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit'][0])) {
					
					foreach ($result['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit'] AS $hit) {
						
						$blastResult[] = $this->hitInfo($hit, $result['BlastOutput_query-len'], isset($request['allowmore']));
						
					}

				} else {
					$hit = $result['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit'];					
					$blastResult[] = $this->hitInfo($hit, $result['BlastOutput_query-len']);					
				}
			}
			return $blastResult;
			
			
		}
	
		function hitInfo($hit, $queryLength, $moreSeq = false)
		{
			$totalScore = 0;
			
			if(isset($hit['Hit_hsps']['Hsp'][0])) {
				foreach($hit['Hit_hsps']['Hsp'] AS $hsp) { 
					$totalScore += $hsp['Hsp_bit-score']; 
				}	
			} else if (!isset($hit['Hit_hsps']['Hsp'][0]) && !empty($hit['Hit_hsps']['Hsp'])) {
				$totalScore = $hit['Hit_hsps']['Hsp']['Hsp_bit-score']; 
			}
			
			$topHsp = (isset($hit['Hit_hsps']['Hsp'][0]))?$hit['Hit_hsps']['Hsp'][0]:$hit['Hit_hsps']['Hsp'];
			
			$maxScore = $topHsp['Hsp_bit-score'];
			
			if($moreSeq) { 
				$itemId = $hit['Hit_accession']; $itemTitle = $hit['Hit_def']; 
			} else { 
				$itemId = explode('_', $hit['Hit_def']); 
				$itemTitle = $itemId[1]; 
				$itemId = $itemId[0]; 
			}
			
			$get_data=DB::table('barcode')->where('barcode_id',$itemId)->get();
			$gene = $get_data[0]['gene'];
			
			$queryCoverage = round((($topHsp['Hsp_align-len']-$topHsp['Hsp_gaps'])/$queryLength)*100, 2).'%';
			$evalue = $topHsp['Hsp_evalue'];
			$ident = sprintf("%d/%d (%01.2f%%)", $topHsp['Hsp_identity'], $topHsp['Hsp_query-to'], round(($topHsp['Hsp_identity']/$topHsp['Hsp_query-to'])*100, 2));

			$alignTop =  str_split($topHsp['Hsp_qseq'], 60);
			$alignMiddle =  str_split($topHsp['Hsp_midline'], 60);
			$alignBottom =  str_split($topHsp['Hsp_hseq'], 60);
			$alignLength = $topHsp['Hsp_align-len'];
			$queryFrom = $topHsp['Hsp_query-from'];
			$queryTo = $topHsp['Hsp_query-to'];
			$hitLength = $hit['Hit_len'];
			$hitFrom = $topHsp['Hsp_hit-from'];
			$hitTo = $topHsp['Hsp_hit-to'];
			$bitScore = $topHsp['Hsp_bit-score'];
			$gaps = $topHsp['Hsp_gaps'];
			$hitSpace = ($hitFrom>$hitTo)?-60:60;
			$totalLine = ceil($alignLength/60);
			
			return $hitInfo = array(
					'maxScore'	=> $maxScore,
					'totalScore'	=> $totalScore,
					'bitScore'	=>	$bitScore,
					'hitLength'	=>	$hitLength,
					'gaps'	=>	$gaps,
					'itemId'	=>	$itemId,
					'gene'		=>	$gene,
					'itemTitle'	=>	$itemTitle,
					'queryCoverage'	=>	$queryCoverage,
					'ident'	=>	$ident,
					'evalue'	=>	$evalue,
					'alignTop'	=>	$alignTop,
					'alignMiddle'	=>	$alignMiddle,
					'alignBottom'	=>	$alignBottom,
					'alignLength'	=>	$alignLength,
					'queryFrom'	=>	$queryFrom,
					'queryTo'	=>	$queryTo,
					'hitFrom'	=>	$hitFrom,
					'hitTo'	=>	$hitTo,
					'hitSpace'	=>	$hitSpace,
					'totalLine'	=>	$totalLine
				);
			
		}
	
	}
	
