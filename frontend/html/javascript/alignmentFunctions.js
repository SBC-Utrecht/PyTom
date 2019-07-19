/**
 * 
 */

function updateAlignmentCookiesFromPage(document){
	//Subtomogram information
	var pixelSize 			= document.getElementById('pixelSize').value;
	var particleDiameter 	= document.getElementById('particleDiameter').value;
	
	setCookie('pytom_pixelSize',pixelSize);
	setCookie('pytom_particleDiameter',particleDiameter);
	
	
	//Sampling settings
	var samplingGlobalRadio = document.getElementById('SamplingGlobal');
	
	var samplingMethod = 'Local';
	
	if(samplingGlobalRadio.checked){
		samplingMethod = 'Global'; 
		var angleList = document.getElementById('AngleList');
		var globalSamplingFile = angleList[angleList.selectedIndex].value;
		setCookie('pytom_angleFile',globalSamplingFile,false);
	}else{
		var angleIncrement 	= document.getElementById('AngleStartIncrement');
		var shellIncrement 	= document.getElementById('AngleShellIncrement');
		var shells 			= document.getElementById('AngleNumberShells');
		
		setCookie('pytom_angleStartIncrement',angleIncrement.value);
		setCookie('pytom_angleShells',shells.value);
		setCookie('pytom_shellIncrement',shellIncrement.value);
	}
	
	setCookie('pytom_samplingMethod',samplingMethod,false);
	
	//filter settings
	var lowestFrequency 	= document.getElementById('LowestFrequency').value;
	var highestFrequency 	= document.getElementById('HighestFrequency').value;
	var smooth				= document.getElementById('FilterSmooth').value;
	
	setCookie('pytom_lowestFrequency',lowestFrequency);
	setCookie('pytom_highestFrequency',highestFrequency);
	setCookie('pytom_filterSmooth',smooth);
	
	//adaptive settings
	var adaptiveRadio = document.getElementById('AdaptiveSamplingOn');
	
	if(adaptiveRadio.checked){
		setCookie('pytom_adaptive','On');
		
		var pytom_resolutionCriterion 	= document.getElementById('CCutoff').value;
		var resolutionOffset 			= document.getElementById('deltaR').value;
		var angleFactor					= document.getElementById('angleFactor').value;
		
		setCookie('pytom_resolutionCriterion',pytom_resolutionCriterion);
		setCookie('pytom_resolutionOffset',resolutionOffset);
		setCookie('pytom_angleFactor',angleFactor);
		
	}else{
		setCookie('pytom_adaptive','Off');
	}
	
	//score
	var scoreType = document.getElementById('ScoreList');
	scoreType = scoreType[scoreType.selectedIndex].value;
	setCookie('pytom_scoreType',scoreType);
	
	//peak prior
	var peakPriorRadius = document.getElementById('peakPriorRadius').value;
	setCookie('pytom_peakPriorRadius',peakPriorRadius);
	
	var peakPriorSmooth = document.getElementById('peakPriorSmooth').value;
	setCookie('pytom_peakPriorSmooth',peakPriorSmooth);
	
	//processing info
	var numberIterations = document.getElementById('numberIterations').value;
	setCookie('pytom_numberIterations',numberIterations);
	
	var binning = document.getElementById('binning').value;
	setCookie('pytom_binning',binning);
}



function checkAlignmentParameters(theOverlay){
	//display overlay
	//request content when clicked, fill destination with content
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){		
		
		var emptyLine 			= "<tr><td align='center' align='middle'><font class='checkTextName'></font></td>"+
								  "<td><font class='checkTextValue'>&nbsp;</font></td>" +
								  "</tr>"
								  
		var particleListDIR		= getCookie('pytom_particleListDIR');
		var particleListXML		= getCookie('pytom_particleListXML');
		var pixelSize			= getCookie('pytom_pixelSize');
		var particleDiameter	= getCookie('pytom_particleDiameter');
		var reference		 	= getCookie('pytom_referenceFile');
		var mask				= getCookie('pytom_maskFile');
		var samplingMethod	    = getCookie('pytom_samplingMethod');
		var angleFile 			= getCookie('pytom_angleFile');
		var angleStart 			= getCookie('pytom_angleStartIncrement');
		var angleShells 		= getCookie('pytom_angleShells');
		var shellIncrement		= getCookie('pytom_shellIncrement');
		var lowestFrequency 	= getCookie('pytom_lowestFrequency');
		var highestFrequency 	= getCookie('pytom_highestFrequency');
		var filterSmooth		= getCookie('pytom_filterSmooth');
		var adaptive			= getCookie('pytom_adaptive');
		var resolutionCriterion	= getCookie('pytom_resolutionCriterion');
		var resolutionOffset	= getCookie('pytom_resolutionOffset');
		var angleFactor			= getCookie('pytom_angleFactor');
		var scoreType			= getCookie('pytom_scoreType');
		var peakPriorRadius		= getCookie('pytom_peakPriorRadius');
		var peakPriorSmooth		= getCookie('pytom_peakPriorSmooth');
		var numberIterations	= getCookie('pytom_numberIterations');
		var binning			= getCookie('pytom_binning');
		var destinationDirectory= getCookie('pytom_destination');
		
		theOverlay.style.display = "";
		
		content = 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Particle List</font></td>"
		
		if(particleListXML != null ){
			"<td><font class='checkTextValue'>&nbsp;" + particleListXML + "</font></td> </tr>";
		}
		else if(particleListDIR != null){
			"<td><font class='checkTextValue'>&nbsp;" + particleListDIR + "</font></td> </tr>";
		}
		else{
			content += "<td><font class='wrongValue'>&nbsp; No particle list selected!</font></td></tr>";
		}
		
					
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Pixel Size</font></td>"+
					"<td><font class='checkTextValue'>&nbsp;"+ pixelSize +"</font></td>" +
					"</tr>"
					
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Partice Diameter</font></td>"+
					"<td><font class='checkTextValue'>&nbsp;"+ particleDiameter +"</font></td>" +
					"</tr>";
		
		content +=  emptyLine;
					
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Reference</font></td>"+
					"<td><font class='checkTextValue'>&nbsp;" + reference + "</font></td>" +
					"</tr>" +
					"<tr><td align='center' align='middle'><font class='checkTextName'>Mask</font></td>"+
					"<td><font class='checkTextValue'>&nbsp;" + mask + "</font></td>" +
					"</tr>"; 
					
	    content +=  emptyLine;
	    
		content +=	"<tr><td align='center' align='middle'><font class='checkTextName'>Sampling Method</font></td>"+
					"<td><font class='checkTextValue'>" + samplingMethod + "</font></td>" +
					"</tr>";
					
		if(samplingMethod == 'Global'){
			
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Sampling List</font></td>"+
						"<td><font class='checkTextValue'>" + angleFile + "</font></td>" +
						"</tr>" 
						
		}else{
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Number Shells</font></td>"+
			"<td><font class='checkTextValue'>" + angleShells + "</font></td>" +
			"</tr>"
			
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Angle Increment</font></td>"+
			"<td><font class='checkTextValue'>" + angleStart + "</font></td>" +
			"</tr>"
			
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Shell Increment</font></td>"+
			"<td><font class='checkTextValue'>" + shellIncrement + "</font></td>" +
			"</tr>"
		}
		
		content +=  emptyLine;
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Lowest Frequency</font></td>"+
		"<td><font class='checkTextValue'>" + lowestFrequency + "</font></td>" +
		"</tr>";
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Highest Frequency</font></td>"+
		"<td><font class='checkTextValue'>" + highestFrequency + "</font></td>" +
		"</tr>";
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Filter Smooth</font></td>"+
		"<td><font class='checkTextValue'>" + filterSmooth + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		
		if(adaptive == "On"){
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Adaptive Sampling </font></td>"+
			"<td><font class='checkTextValue'>On</font></td>" +
			"</tr>";
			
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Resolution Criterion </font></td>"+
			"<td><font class='checkTextValue'>"+ resolutionCriterion +"</font></td>" +
			"</tr>";
			
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Resolution Scale </font></td>"+
			"<td><font class='checkTextValue'>"+ resolutionOffset +"</font></td>" +
			"</tr>";
			
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Angle Factor </font></td>"+
			"<td><font class='checkTextValue'>"+ angleFactor +"</font></td>" +
			"</tr>";
			
		}else{
			content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Adaptive Sampling </font></td>"+
			"<td><font class='checkTextValue'>Off</font></td>" +
			"</tr>";
		}
		
		content +=  emptyLine;
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Score </font></td>"+
		"<td><font class='checkTextValue'>"+ scoreType +"</font></td>" +
		"</tr>";
		content +=  emptyLine;
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Peak Prior Radius </font></td>"+
		"<td><font class='checkTextValue'>"+ peakPriorRadius +"</font></td>" +
		"</tr>";
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Peak Prior Smooth </font></td>"+
		"<td><font class='checkTextValue'>"+ peakPriorSmooth +"</font></td>" +
		"</tr>";
		content +=  emptyLine;
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Number iterations </font></td>"+
		"<td><font class='checkTextValue'>"+ numberIterations +"</font></td>" +
		"</tr>";
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Binning </font></td>"+
		"<td><font class='checkTextValue'>"+ binning +"</font></td>" +
		"</tr>";
		content +=  emptyLine;
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Result dir </font></td>"+
		"<td><font class='checkTextValue'>"+ destinationDirectory +"</font></td>" +
		"</tr>";
		
		content += 	"<tr><td></td><td><input type='button' value='Cancel' onClick='return hideOverlay()'/></td></tr></table></center>";
							   
		theOverlay.innerHTML = content; 
	}else{
		hideOverlay();
	}
	return false;
}


function createAlignmentCookies(){
	
	setCookie('pytom_particleListXML','',false);
	setCookie('pytom_referenceFile','',false);
	setCookie('pytom_maskFile','',false);
	setCookie('pytom_samplingMethod','',false);
	setCookie('pytom_angleFile','',false);
	setCookie('pytom_angleStartIncrement','',false);
	setCookie('pytom_angleShells','',false);
	setCookie('pytom_lowestFrequency','',false);
	setCookie('pytom_highestFrequency','',false);
	setCookie('pytom_smoothFrequency','',false);
	setCookie('pytom_adaptive','',false);
	setCookie('pytom_resolutionCriterion','',false);
	setCookie('pytom_resolutionOffset','',false);
	setCookie('pytom_angleFactor','',false);
	setCookie('pytom_pixelSize','',false);
	setCookie('pytom_particleDiameter','',false);
	setCookie('pytom_scoreType','',false);
	setCookie('pytom_peakPriorRadius','',false);
	setCookie('pytom_peakPriorRadius','',false);
	setCookie('pytom_numberIterations','',false);
	setCookie('pytom_binning','',false);
	setCookie('pytom_destination','',false);
}

function eraseAlignmentCookies(){
	
	setCookie('pytom_particleListXML','',true);
	setCookie('pytom_referenceFile','',true);
	setCookie('pytom_maskFile','',true);
	setCookie('pytom_samplingMethod','',true);
	setCookie('pytom_angleFile','',true);
	setCookie('pytom_angleStartIncrement','',true);
	setCookie('pytom_angleShells','',true);
	setCookie('pytom_lowestFrequency','',true);
	setCookie('pytom_highestFrequency','',true);
	setCookie('pytom_smoothFrequency','',true);
	setCookie('pytom_adaptive','',true);
	setCookie('pytom_resolutionCriterion','',true);
	setCookie('pytom_resolutionOffset','',true);
	setCookie('pytom_angleFactor','',true);
	setCookie('pytom_pixelSize','',true);
	setCookie('pytom_particleDiameter','',true);
	setCookie('pytom_scoreType','',true);
	setCookie('pytom_peakPriorRadius','',true);
	setCookie('pytom_peakPriorRadius','',true);
	setCookie('pytom_numberIterations','',true);
	setCookie('pytom_binning','',true);
	setCookie('pytom_destination','',true);
}

function createAlignmentJobFromParameters(){
	
	var particleListXML		= getCookie('pytom_particleListXML');
	var particleListDIR		= getCookie('pytom_particleListDIR');
	var pixelSize			= getCookie('pytom_pixelSize');
	var particleDiameter	= getCookie('pytom_particleDiameter');
	var reference		 	= getCookie('pytom_referenceFile');
	var mask				= getCookie('pytom_maskFile');
	var samplingMethod	    = getCookie('pytom_samplingMethod');
	var angleFile 			= getCookie('pytom_angleFile');
	var angleStart 			= getCookie('pytom_angleStartIncrement');
	var angleShells 		= getCookie('pytom_angleShells');
	var shellIncrement		= getCookie('pytom_shellIncrement');
	var lowestFrequency 	= getCookie('pytom_lowestFrequency');
	var highestFrequency 	= getCookie('pytom_highestFrequency');
	var filterSmooth		= getCookie('pytom_filterSmooth');
	var adaptive			= getCookie('pytom_adaptive');
	var resolutionCriterion	= getCookie('pytom_resolutionCriterion');
	var resolutionOffset	= getCookie('pytom_resolutionOffset');
	var angleFactor			= getCookie('pytom_angleFactor');
	var scoreType			= getCookie('pytom_scoreType');
	var peakPriorRadius		= getCookie('pytom_peakPriorRadius');
	var peakPriorSmooth		= getCookie('pytom_peakPriorSmooth');
	var numberIterations	= getCookie('pytom_numberIterations');
	var binning			= getCookie('pytom_binning');
	var destinationDirectory= getCookie('pytom_destination');
	
	var requestString = 'createAlignmentJob.py?';

	if(particleListXML != null){
		requestString += 'plXML=' + particleListXML + '&';
	}
	else if(particleListDIR != null){
		requestString += 'plDIR=' + particleListDIR + '&';
	}
	
	//pixel size and particle diameter
	requestString += 'pixSize=' + pixelSize + '&';
	requestString += 'partDia=' + particleDiameter+ '&';
	requestString += 'ref=' + reference + '&';
	requestString += 'mask=' + mask + '&';
	
	if(samplingMethod == 'Global'){
		requestString += 'sampling=GLOBAL&';
		requestString += 'angFile=' + angleFile+ '&';
	}else{
		requestString += 'sampling=LOCAL&';
		requestString += 'angStart=' + angleStart+ '&';
		requestString += 'angShells=' + angleShells+ '&';
		requestString += 'angInc=' + shellIncrement+ '&';
	}
	
	requestString += 'lowestF=' + lowestFrequency + '&';
	requestString += 'highestF=' + highestFrequency + '&';
	requestString += 'filtSm=' + filterSmooth + '&';
	
	if(adaptive == 'On'){
		requestString += 'adapt=ON&';
		requestString += 'adResC=' + resolutionCriterion + '&';
		requestString += 'adResOf=' + resolutionOffset + '&';
		requestString += 'angFac=' + angleFactor + '&';
	}else{
		requestString += 'adapt=OFF&';
	}
	
	requestString += 'score=' + scoreType + '&';
	
	requestString += 'pkPriRad=' + peakPriorRadius + '&';
	requestString += 'pkSmooth=' + peakPriorSmooth + '&';
	
	requestString += 'iter=' + numberIterations + '&';
	requestString += 'binning=' + binning + '&';
	requestString += 'dest=' + destinationDirectory;

	//display some additional boxes like - set filename of job, ...
	jobFile = prompt("Choose a Jobfile name. Provide full path! End with .xml");
	
	if(jobFile.length > 0){
		requestString += '&jobFile=' + jobFile;
	} 
	
	//send request to server
	request(requestString,true);
	
	
	//display success or fail
	
}



function resetAlignmentValues(){
	
	//Subtomogram information
	document.getElementById('pixelSize').value = '';
	document.getElementById('particleDiameter').value = '';
	
	//Sampling settings
	var samplingGlobalRadio = document.getElementById('SamplingGlobal');
	
	if(samplingGlobalRadio.checked){
		samplingGlobalRadio.checked = false; 
	}else{

	}
	
	//filter settings
	document.getElementById('LowestFrequency').value = '0';
	document.getElementById('HighestFrequency').value = '0';
	document.getElementById('FilterSmooth').value = '0';
	
	//adaptive settings
	var adaptiveRadio = document.getElementById('AdaptiveSamplingOn');
	
	if(adaptiveRadio.checked){
		adaptiveRadio.checked = false;
		
	}else{
	
	}
	
	document.getElementById('CCutoff').value = '0.5';
	document.getElementById('deltaR').value = 0.1;
	document.getElementById('angleFactor').value = 0.5;
	
	//score
	var scoreType = document.getElementById('ScoreList');
	scoreType[scoreType.selectedIndex].value = 'XCF';
	
	//peak prior
	document.getElementById('peakPriorRadius').value = '-1';
	document.getElementById('peakPriorSmooth').value = '-1';
	
	//processing info
	document.getElementById('numberIterations').value = '10';
	document.getElementById('binning').value = '1';
	
}


function loadAlignmentJob(){
	overlayWaitMessageRunFunction(document.getElementById('overlayBox'),'Select <b>existing Alignment Job</b><br/>','loadAlignmentJob.py','File',alert)
}

function parseAlignmentJobResponse(){
	
}


