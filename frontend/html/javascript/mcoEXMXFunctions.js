/**
 * 
 */

function updateMCOEXMXCookiesFromPage(document){
	var lowestFrequency 	= document.getElementById('LowestFrequency').value;
	var highestFrequency 	= document.getElementById('HighestFrequency').value;
	var smooth				= document.getElementById('FilterSmooth').value;
	
	setCookie('pytom_mcoEXMX_lowestFrequency',lowestFrequency);
	setCookie('pytom_mcoEXMX_highestFrequency',highestFrequency);
	setCookie('pytom_mcoEXMX_filterSmooth',smooth);
	
	var scoreType 		= document.getElementById('ScoreList');
	scoreType = scoreType[scoreType.selectedIndex].value;
	setCookie('pytom_mcoEXMX_scoreType',scoreType);
	
	//processing info
	var numberIterations= document.getElementById('numberIterations').value;
	setCookie('pytom_mcoEXMX_numberIterations',numberIterations);
	
	var binning 		= document.getElementById('binning').value;
	setCookie('pytom_mcoEXMX_binning',binning);
	
	var convergence 	= document.getElementById('convergence').value;
	setCookie('pytom_mcoEXMX_convergence',convergence);
	
	var numberClasses 	= document.getElementById('numberClasses').value;
	setCookie('pytom_mcoEXMX_numberClasses',numberClasses);
	
	var wedgeAngle1 	= document.getElementById('Wedge1').value;
	setCookie('pytom_mcoEXMX_wedgeAngle1',wedgeAngle1);
	var wedgeAngle2 	= document.getElementById('Wedge2').value;
	setCookie('pytom_mcoEXMX_wedgeAngle2',wedgeAngle2);
	
	var pixelSize 	= document.getElementById('pixelSize').value;
	setCookie('pytom_mcoEXMX_pixelSize',pixelSize);
	var particleDiameter 	= document.getElementById('particleDiameter').value;
	setCookie('pytom_mcoEXMX_particleDiameter',particleDiameter);
	
}

function createMCOEXMXCookies(){
	setCookie('pytom_mcoEXMX_particleListXML','',false);
	setCookie('pytom_mcoEXMX_particleListDIR','',false);
	setCookie('pytom_mcoEXMX_maskFile','',false);
	setCookie('pytom_mcoEXMX_lowestFrequency','',false);
	setCookie('pytom_mcoEXMX_highestFrequency','',false);
	setCookie('pytom_mcoEXMX_filterSmooth','',false);
	setCookie('pytom_mcoEXMX_scoreType','',false);
	setCookie('pytom_mcoEXMX_numberIterations','',false);
	setCookie('pytom_mcoEXMX_numberClasses','',false);
	setCookie('pytom_mcoEXMX_convergence','',false);
	setCookie('pytom_mcoEXMX_binning','',false);
	setCookie('pytom_mcoEXMX_destination','',false);
	setCookie('pytom_mcoEXMX_pixelSize','',false);
	setCookie('pytom_mcoEXMX_particleDiameter','',false);
	setCookie('pytom_mcoEXMX_wedgeAngle1','',false);
	setCookie('pytom_mcoEXMX_wedgeAngle2','',false);
	
}

function eraseMCOEXMXCookies(){
	setCookie('pytom_mcoEXMX_particleListXML','',true);
	setCookie('pytom_mcoEXMX_particleListDIR','',true);
	setCookie('pytom_mcoEXMX_maskFile','',true);
	setCookie('pytom_mcoEXMX_lowestFrequency','',true);
	setCookie('pytom_mcoEXMX_highestFrequency','',true);
	setCookie('pytom_mcoEXMX_filterSmooth','',true);
	setCookie('pytom_mcoEXMX_scoreType','',true);
	setCookie('pytom_mcoEXMX_numberIterations','',true);
	setCookie('pytom_mcoEXMX_numberClasses','',true);
	setCookie('pytom_mcoEXMX_destination','',true);
	setCookie('pytom_mcoEXMX_convergence','',true);
	setCookie('pytom_mcoEXMX_binning','',true);
	setCookie('pytom_mcoEXMX_pixelSize','',true);
	setCookie('pytom_mcoEXMX_particleDiameter','',true);
	setCookie('pytom_mcoEXMX_wedgeAngle1','',true);
	setCookie('pytom_mcoEXMX_wedgeAngle2','',true);
	
}

function resetMCOEXMXValues(document){
	//reset bandpass
	document.getElementById('LowestFrequency').value = '0';
	document.getElementById('HighestFrequency').value = '1';
	document.getElementById('FilterSmooth').value = '3';
	
	//reset particle info
	document.getElementById('pixelSize').value = '0';
	document.getElementById('particleDiameter').value = '0';
	
	//reset iteration settings
	document.getElementById('numberIterations').value = '';
	document.getElementById('binning').value = '';
	document.getElementById('convergence').value = '';
	
	//reset wedge settings
	document.getElementById('Wedge1').value = '';
	document.getElementById('Wedge2').value = '';
	
	document.getElementById('numberClasses').value = '2';
}


function checkMCOEXMXParameters(theOverlay){
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		
		theOverlay.style.display = "";
		
		var particleListDIR		= getCookie('pytom_mcoEXMX_particleListDIR');
		var particleListXML		= getCookie('pytom_mcoEXMX_particleListXML');
		var mask				= getCookie('pytom_mcoEXMX_maskFile');
		var lowestFrequency 	= getCookie('pytom_mcoEXMX_lowestFrequency');
		var highestFrequency 	= getCookie('pytom_mcoEXMX_highestFrequency');
		var filterSmooth		= getCookie('pytom_mcoEXMX_filterSmooth');
		var numberIterations	= getCookie('pytom_mcoEXMX_numberIterations');
		var numberClasses		= getCookie('pytom_mcoEXMX_numberClasses');
		var binning			= getCookie('pytom_mcoEXMX_binning');
		var convergence 		= getCookie('pytom_mcoEXMX_convergence');
		var destinationDirectory= getCookie('pytom_mcoEXMX_destination');
		var pixelSize			= getCookie('pytom_mcoEXMX_pixelSize');
		var particleDiameter	= getCookie('pytom_mcoEXMX_particleDiameter');
		var wedgeAngle1 		= getCookie('pytom_mcoEXMX_wedgeAngle1');
		var wedgeAngle2 		= getCookie('pytom_mcoEXMX_wedgeAngle2');
		var scoreType 			= getCookie('pytom_mcoEXMX_scoreType');
		
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
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Number of Iterations</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + numberIterations + "</font></td>" +
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
		
		content += 	"<tr><td></td><td><input type='button' value='Cancel' onClick='return hideOverlay()'/></td></tr></table></center>";
		   
		theOverlay.innerHTML = content; 
	}else{
		hideOverlay();
	}
	return false;
	
}



function creatMCOEXMXJobFromParameters(theOverlay){

	var particleListDIR		= getCookie('pytom_mcoEXMX_particleListDIR');
	var particleListXML		= getCookie('pytom_mcoEXMX_particleListXML');
	var mask				= getCookie('pytom_mcoEXMX_maskFile');
	var lowestFrequency 	= getCookie('pytom_mcoEXMX_lowestFrequency');
	var highestFrequency 	= getCookie('pytom_mcoEXMX_highestFrequency');
	var filterSmooth		= getCookie('pytom_mcoEXMX_filterSmooth');
	var numberIterations	= getCookie('pytom_mcoEXMX_numberIterations');
	var numberClasses		= getCookie('pytom_mcoEXMX_numberClasses');
	var binning			= getCookie('pytom_mcoEXMX_binning');
	var convergence 		= getCookie('pytom_mcoEXMX_convergence');
	var destinationDirectory= getCookie('pytom_mcoEXMX_destination');
	var pixelSize			= getCookie('pytom_mcoEXMX_pixelSize');
	var particleDiameter	= getCookie('pytom_mcoEXMX_particleDiameter');
	var wedgeAngle1 		= getCookie('pytom_mcoEXMX_wedgeAngle1');
	var wedgeAngle2 		= getCookie('pytom_mcoEXMX_wedgeAngle2');
	var scoreType 			= getCookie('pytom_mcoEXMX_scoreType');
	
	
	var requestString = 'createMCOEXMXJob.py?';
	
	if(particleListXML != null){
		requestString += 'plXML=' + particleListXML + '&';
	}
	else if(particleListDIR != null){
		requestString += 'plDIR=' + particleListDIR + '&';
	}
	
	//pixel size and particle diameter
	requestString += 'pixSize='		+ pixelSize + '&';
	requestString += 'partDia=' 	+ particleDiameter+ '&';
	requestString += 'mask=' 		+ mask + '&';
	
	requestString += 'lowestF=' 	+ lowestFrequency + '&';
	requestString += 'highestF=' 	+ highestFrequency + '&';
	requestString += 'filtSm=' 		+ filterSmooth + '&';
	
	requestString += 'score=' 	+ scoreType + '&';
	
	requestString += 'iter=' 	+ numberIterations + '&';
	requestString += 'classes=' + numberClasses + '&';
	requestString += 'binning=' 	+ binning + '&';
	requestString += 'dest=' 	+ destinationDirectory + '&';
	requestString += 'conv=' 	+ convergence + '&';
	
	requestString += 'wa1=' + wedgeAngle1 + '&';
	requestString += 'wa2=' + wedgeAngle2 + '&';
	
	
	//display some additional boxes like - set filename of job, ...
	jobFile = prompt("Choose a Jobfile name. Provide full path! End with XML!");
	
	if(jobFile.length > 0){
		requestString += 'jobFile=' + jobFile;
	} 
	
	//send request to server
	request(requestString,true);
}