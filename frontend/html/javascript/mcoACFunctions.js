/**
 * 
 */

function updateMCOACCookiesFromPage(){
	var lowestFrequency 	= document.getElementById('LowestFrequency').value;
	var highestFrequency 	= document.getElementById('HighestFrequency').value;
	var smooth				= document.getElementById('FilterSmooth').value;
	
	setCookie('pytom_mcoAC_lowestFrequency',lowestFrequency);
	setCookie('pytom_mcoAC_highestFrequency',highestFrequency);
	setCookie('pytom_mcoAC_filterSmooth',smooth);
	
	var scoreType 		= document.getElementById('ScoreList');
	scoreType = scoreType[scoreType.selectedIndex].value;
	setCookie('pytom_mcoAC_scoreType',scoreType);
	
	var binning 		= document.getElementById('binning').value;
	setCookie('pytom_mcoAC_binning',binning);
	
	var convergence 	= document.getElementById('convergence').value;
	setCookie('pytom_mcoAC_convergence',convergence);
	
	var numberClasses 	= document.getElementById('numberClasses').value;
	setCookie('pytom_mcoAC_numberClasses',numberClasses);
	
	var wedgeAngle1 	= document.getElementById('Wedge1').value;
	setCookie('pytom_mcoAC_wedgeAngle1',wedgeAngle1);
	var wedgeAngle2 	= document.getElementById('Wedge2').value;
	setCookie('pytom_mcoAC_wedgeAngle2',wedgeAngle2);
	
	var pixelSize 			= document.getElementById('pixelSize').value;
	setCookie('pytom_mcoAC_pixelSize',pixelSize);
	
	var particleDiameter 	= document.getElementById('particleDiameter').value;
	setCookie('pytom_mcoAC_particleDiameter',particleDiameter);
	
	var numberRefinements	= document.getElementById('numberRefinements').value;
	setCookie('pytom_mcoAC_numberRefinements',numberRefinements);
	
	var temperatureType 	= document.getElementById('temperatureType');
	setCookie('pytom_mcoAC_temperature',temperatureType[temperatureType.selectedIndex].value);
	
	var classificationCriterion 	= document.getElementById('classificationCriterion');
	setCookie('pytom_mcoAC_criterion',classificationCriterion[classificationCriterion.selectedIndex].value);
	
	var startTemperature	= document.getElementById('startTemperature').value;
	setCookie('pytom_mcoAC_startTemperature',startTemperature);
	
	var annealingStep	= document.getElementById('annealingStep').value;
	setCookie('pytom_mcoAC_annealingStep',annealingStep);
	
}

function createMCOACCookies(){
	setCookie('pytom_mcoAC_particleListXML','',false);
	setCookie('pytom_mcoAC_particleListDIR','',false);
	setCookie('pytom_mcoAC_maskFile','',false);
	setCookie('pytom_mcoAC_lowestFrequency','',false);
	setCookie('pytom_mcoAC_highestFrequency','',false);
	setCookie('pytom_mcoAC_filterSmooth','',false);
	setCookie('pytom_mcoAC_scoreType','',false);
	setCookie('pytom_mcoAC_numberClasses','',false);
	setCookie('pytom_mcoAC_convergence','',false);
	setCookie('pytom_mcoAC_binning','',false);
	setCookie('pytom_mcoAC_destination','',false);
	setCookie('pytom_mcoAC_pixelSize','',false);
	setCookie('pytom_mcoAC_particleDiameter','',false);
	setCookie('pytom_mcoAC_wedgeAngle1','',false);
	setCookie('pytom_mcoAC_wedgeAngle2','',false);
	setCookie('pytom_mcoAC_temperature','',false);
	setCookie('pytom_mcoAC_criterion','',false);
	setCookie('pytom_mcoAC_numberRefinements','',false);
	setCookie('pytom_mcoAC_startTemperature','',false);
	setCookie('pytom_mcoAC_annealingStep','',false);
}

function eraseMCOACCookies(){
	setCookie('pytom_mcoAC_particleListXML','',true);
	setCookie('pytom_mcoAC_particleListDIR','',true);
	setCookie('pytom_mcoAC_maskFile','',true);
	setCookie('pytom_mcoAC_lowestFrequency','',true);
	setCookie('pytom_mcoAC_highestFrequency','',true);
	setCookie('pytom_mcoAC_filterSmooth','',true);
	setCookie('pytom_mcoAC_scoreType','',true);
	setCookie('pytom_mcoAC_numberClasses','',true);
	setCookie('pytom_mcoAC_destination','',true);
	setCookie('pytom_mcoAC_convergence','',true);
	setCookie('pytom_mcoAC_binning','',true);
	setCookie('pytom_mcoAC_pixelSize','',true);
	setCookie('pytom_mcoAC_particleDiameter','',true);
	setCookie('pytom_mcoAC_wedgeAngle1','',true);
	setCookie('pytom_mcoAC_wedgeAngle2','',true);
	setCookie('pytom_mcoAC_temperature','',true);
	setCookie('pytom_mcoAC_criterion','',true);
	setCookie('pytom_mcoAC_numberRefinements','',true);
	setCookie('pytom_mcoAC_startTemperature','',true);
	setCookie('pytom_mcoAC_annealingStep','',true);
}

function resetMCOACValues(document){
	//reset bandpass
	document.getElementById('LowestFrequency').value = '0';
	document.getElementById('HighestFrequency').value = '1';
	document.getElementById('FilterSmooth').value = '3';
	
	//reset particle info
	document.getElementById('pixelSize').value = '0';
	document.getElementById('particleDiameter').value = '0';
	
	//reset iteration settings
	document.getElementById('binning').value = '';
	document.getElementById('convergence').value = '';
	document.getElementById('numberClasses').value = '2';
	
	//reset wedge settings
	document.getElementById('Wedge1').value = '';
	document.getElementById('Wedge2').value = '';
	
	document.getElementById('numberRefinements').value 	= '';
	document.getElementById('startTemperature').value 	= '1';
	document.getElementById('annealingStep').value 	= '0.1';
}


function checkMCOACParameters(theOverlay){
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		
		theOverlay.style.display = "";
		
		var particleListDIR		= getCookie('pytom_mcoAC_particleListDIR');
		var particleListXML		= getCookie('pytom_mcoAC_particleListXML');
		var mask				= getCookie('pytom_mcoAC_maskFile');
		var lowestFrequency 	= getCookie('pytom_mcoAC_lowestFrequency');
		var highestFrequency 	= getCookie('pytom_mcoAC_highestFrequency');
		var filterSmooth		= getCookie('pytom_mcoAC_filterSmooth');
		var numberClasses		= getCookie('pytom_mcoAC_numberClasses');
		var binning			= getCookie('pytom_mcoAC_binning');
		var convergence 		= getCookie('pytom_mcoAC_convergence');
		var destinationDirectory= getCookie('pytom_mcoAC_destination');
		var pixelSize			= getCookie('pytom_mcoAC_pixelSize');
		var particleDiameter	= getCookie('pytom_mcoAC_particleDiameter');
		var wedgeAngle1 		= getCookie('pytom_mcoAC_wedgeAngle1');
		var wedgeAngle2 		= getCookie('pytom_mcoAC_wedgeAngle2');
		var scoreType 			= getCookie('pytom_mcoAC_scoreType');
		var temperature 		= getCookie('pytom_mcoAC_temperature');
		var criterion 			= getCookie('pytom_mcoAC_criterion');
		var numberRefinements	= getCookie('pytom_mcoAC_numberRefinements');
		var startTemperature	= getCookie('pytom_mcoAC_startTemperature');
		var annealingStep		= getCookie('pytom_mcoAC_annealingStep');
		var particleList = '';
		
		if(particleListXML != null){
			particleList = particleListXML;
		}
		else if(particleListDIR != null){
			particleList = particleListDIR;
		}
		
		if(particleList == ''){
			particleList = 'Not set'
		}
		
		var emptyLine 			= "<tr><td align='center' align='middle'></td><td></td></tr>"
		
		var content = 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Particle List</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + particleList + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Mask</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + mask + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Score </font></td>"+
		"<td><font class='checkTextValue'>" + scoreType +"</font></td>" +
		"</tr>";
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
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Wedge Angle 1</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + wedgeAngle1 + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Wedge Angle 2</font></td>"+
		"<td><font class='checkTextValue'>" + wedgeAngle2 +"</font></td>" +
		"</tr>";
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Number of classes</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + numberClasses + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content += 	"<tr><td align='center' align='middle'><font class='checkTextName'>Binning</font></td>"+
		"<td><font class='checkTextValue'>" + binning +"</font></td>" +
		"</tr>";
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Convergence criterion</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + convergence + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Destination directory</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + destinationDirectory + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Pixel size</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + pixelSize + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Particle diameter</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + particleDiameter + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Temperature type</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + temperature + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;

		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Start temperature</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + startTemperature + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Annealing temperature step</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + annealingStep + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Annealing criterion</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + criterion + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Number local refinements</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + numberRefinements + "</font></td>" +
		"</tr>"; 
		
		content +=  emptyLine;
		
		
		content += 	"<tr><td></td><td><input type='button' value='Cancel' onClick='return hideOverlay()'/></td></tr></table></center>";
		   
		theOverlay.innerHTML = content; 
	}else{
		hideOverlay();
	}
	return false;
	
}



function creatMCOACJobFromParameters(){
	updateMCOACCookiesFromPage();
	var particleListDIR		= getCookie('pytom_mcoAC_particleListDIR');
	var particleListXML		= getCookie('pytom_mcoAC_particleListXML');
	var mask				= getCookie('pytom_mcoAC_maskFile');
	var lowestFrequency 	= getCookie('pytom_mcoAC_lowestFrequency');
	var highestFrequency 	= getCookie('pytom_mcoAC_highestFrequency');
	var filterSmooth		= getCookie('pytom_mcoAC_filterSmooth');
	var numberClasses		= getCookie('pytom_mcoAC_numberClasses');
	var binning			= getCookie('pytom_mcoAC_binning');
	var convergence 		= getCookie('pytom_mcoAC_convergence');
	var destinationDirectory= getCookie('pytom_mcoAC_destination');
	var pixelSize			= getCookie('pytom_mcoAC_pixelSize');
	var particleDiameter	= getCookie('pytom_mcoAC_particleDiameter');
	var wedgeAngle1 		= getCookie('pytom_mcoAC_wedgeAngle1');
	var wedgeAngle2 		= getCookie('pytom_mcoAC_wedgeAngle2');
	var scoreType 			= getCookie('pytom_mcoAC_scoreType');
	var temperature 		= getCookie('pytom_mcoAC_temperature');
	var criterion 			= getCookie('pytom_mcoAC_criterion');
	var numberRefinements	= getCookie('pytom_mcoAC_numberRefinements');
	var startTemperature	= getCookie('pytom_mcoAC_startTemperature');
	var annealingStep		= getCookie('pytom_mcoAC_annealingStep');
	
	var requestString = 'createMCOACJob.py?';
	
	if(particleListXML != null){
		requestString += 'plXML=' + particleListXML + '&';
	}
	else if(particleListDIR != null){
		requestString += 'plDIR=' + particleListDIR + '&';
	}
	
	//pixel size and particle diameter
	requestString += 'pixSize='	+ pixelSize + '&';
	requestString += 'partDia=' 	+ particleDiameter+ '&';
	requestString += 'mask=' 		+ mask + '&';
	
	requestString += 'lowestF=' 	+ lowestFrequency + '&';
	requestString += 'highestF=' 	+ highestFrequency + '&';
	requestString += 'filtSm=' 	+ filterSmooth + '&';
	
	requestString += 'score=' + scoreType + '&';
	
	requestString += 'classes=' + numberClasses + '&';
	requestString += 'binning=' + binning + '&';
	requestString += 'dest=' 	+ destinationDirectory + '&';
	requestString += 'conv=' 	+ convergence + '&';
	
	requestString += 'wa1=' + wedgeAngle1 + '&';
	requestString += 'wa2=' + wedgeAngle2 + '&';
	
	requestString += 'temp=' 		+ temperature		+ '&';
	requestString += 'stemp=' 	+ startTemperature	+ '&';
	requestString += 'astep='		+ annealingStep		+ '&';
	requestString += 'crit=' 		+ criterion			+ '&';
	requestString += 'refin=' 	+ numberRefinements	+ '&';
	
	//display some additional boxes like - set filename of job, ...
	jobFile = prompt("Choose a Jobfile name. Provide full path! End with XML!");
	
	if(jobFile.length > 0){
		requestString += 'jobFile=' + jobFile;
	} 
	
	//send request to server
	request(requestString,true);
}