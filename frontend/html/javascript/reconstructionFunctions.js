/**
 * 
 */
function updateReconstructionCookiesFromPage(){
	
	var xCube = document.getElementById('xCube').value;
	setCookie('pytom_reconstruction_xCube',xCube,false);
	
	var yCube = document.getElementById('yCube').value;
	setCookie('pytom_reconstruction_yCube',yCube,false);
	
	var zCube = document.getElementById('zCube').value;
	setCookie('pytom_reconstruction_zCube',zCube,false);
	
	var binBefore = document.getElementById('binBefore').value;
	setCookie('pytom_reconstruction_binBefore',binBefore,false);
	
	var binAfter = document.getElementById('binAfter').value;
	setCookie('pytom_reconstruction_binAfter',binAfter,false);
	
	try{
		//this function is used to update tomogram and subtomo reconstruction cookies
		//hence, these elements are not available in both html files, thats the reason for the try catch
		var tomogramName = document.getElementById('TomogramFileName').value;
		setCookie('pytom_reconstruction_tomogramName',tomogramName,false);
	}
	catch(e){
	}
	
	try{		
		var xCutPos = document.getElementById('xCutPos').value;
		setCookie('pytom_reconstruction_xCutPos',xCutPos,false);
	}
	catch(e){
	}
	
	try{
		var yCutPos = document.getElementById('yCutPos').value;
		setCookie('pytom_reconstruction_yCutPos',yCutPos,false);
	}
	catch(e){
	}
	
	try{
		var zCutPos = document.getElementById('zCutPos').value;
		setCookie('pytom_reconstruction_zCutPos',zCutPos,false);
	}
	catch(e){
	}
	
	try{
		var tomoBin = document.getElementById('tomoBin').value;
		setCookie('pytom_reconstruction_tomoBin',tomoBin,false);
	}
	catch(e){
	}
}

function eraseReconstructionCookiesFromPage(){
	
	setCookie('pytom_reconstruction_projectionListDIR','',true);
	setCookie('pytom_reconstruction_projectionListXML','',true);
	
	setCookie('pytom_reconstruction_particleListXML','',true);
	
	setCookie('pytom_reconstruction_xCube','',true);
	setCookie('pytom_reconstruction_yCube','',true);
	setCookie('pytom_reconstruction_zCube','',true);
	
	setCookie('pytom_reconstruction_binBefore','',true);
	setCookie('pytom_reconstruction_binAfter','',true);
	
}

function checkReconstructionParameters(theOverlay){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		var emptyLine 			= "<tr><td align='center' align='middle'><font class='checkTextName'></font></td>"+
		  "<td><font class='checkTextValue'>&nbsp;</font></td>" +
		  "</tr>"
		  
		updateReconstructionCookiesFromPage();
		
		var projectionListDIR		= getCookie('pytom_reconstruction_projectionListDIR');
		var projectionListXML		= getCookie('pytom_reconstruction_projectionListXML');
		var particleListXML			= getCookie('pytom_reconstruction_particleListXML');
		var xCube 					= getCookie('pytom_reconstruction_xCube');
		var yCube 					= getCookie('pytom_reconstruction_yCube');
		var zCube 					= getCookie('pytom_reconstruction_zCube');
		var binBefore				= getCookie('pytom_reconstruction_binBefore');
		var binAfter				= getCookie('pytom_reconstruction_binAfter');
		var xCutPos 				= getCookie('pytom_reconstruction_xCutPos');
		var yCutPos 				= getCookie('pytom_reconstruction_yCutPos');
		var zCutPos 				= getCookie('pytom_reconstruction_zCutPos');
		var tomoBin 				= getCookie('pytom_reconstruction_tomoBin');
		
		theOverlay.style.display = "";
		
		content = 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Particle List</font></td>"
		
		if(particleListXML != null){
			content +=  "<td><font class='checkTextValue'>&nbsp;" + particleListXML + "</font></td>" + "</tr>";
		}
		else{
			content += "<td><font class='wrongValue'>&nbsp; No particle list selected!</font></td>" + "</tr>";
		}
		
		content +=  emptyLine;
		
		content += "<tr><td align='center' align='middle'><font class='checkTextName'>Projection List</font></td>"
		if(projectionListXML != null){
			content += "<td><font class='checkTextValue'>&nbsp;" + projectionListXML + "</font></td>" +	"</tr>";
		}
		else if(projectionListDIR != null){
			content += "<td><font class='checkTextValue'>&nbsp;" + projectionListDIR + "</font></td>" + "</tr>";
		}
		else{
			content += "<td><font class='wrongValue'>&nbsp; No projection list selected!</font></td>" +	"</tr>";
		}
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Cube X</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + xCube + "</font></td>" + "</tr>";
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Cube Y</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + yCube + "</font></td>" + "</tr>";
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Cube Z</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + zCube + "</font></td>" + "</tr>";
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Bin before</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + binBefore + "</font></td>" + "</tr>";
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Bin after</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + binAfter + "</font></td>" + "</tr>";
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Cut X</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + xCutPos + "</font></td>" + "</tr>";
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Cut Y</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + yCutPos + "</font></td>" + "</tr>";
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Cut Z</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + zCutPos + "</font></td>" + "</tr>";
		
		content +=  emptyLine;
		
		content +=  "<tr><td align='center' align='middle'><font class='checkTextName'>Tomogram binning</font></td>" +
					"<td><font class='checkTextValue'>&nbsp;" + tomoBin + "</font></td>" + "</tr>";
		
		content += 	"<tr><td></td><td><input type='button' value='Cancel' onClick='return hideOverlay()'/></td></tr></table></center>";
		
		
		theOverlay.innerHTML = content; 
	}else{
		hideOverlay();
	}
	
	return false;
	
}

function createReconstructionJobFromParameters(document){
	
	updateReconstructionCookiesFromPage();
	
	var projectionListDIR		= getCookie('pytom_reconstruction_projectionListDIR');
	var projectionListXML		= getCookie('pytom_reconstruction_projectionListXML');
	var particleListXML			= getCookie('pytom_reconstruction_particleListXML');
	var xCube 					= getCookie('pytom_reconstruction_xCube');
	var yCube 					= getCookie('pytom_reconstruction_yCube');
	var zCube 					= getCookie('pytom_reconstruction_zCube');
	var Before				= getCookie('pytom_reconstruction_Before');
	var After				= getCookie('pytom_reconstruction_After');
	var xCutPos 				= getCookie('pytom_reconstruction_xCutPos');
	var yCutPos 				= getCookie('pytom_reconstruction_yCutPos');
	var zCutPos 				= getCookie('pytom_reconstruction_zCutPos');
	var tomoBin 				= getCookie('pytom_reconstruction_tomoBin');
	
	requestString = 'createReconstructionJob.py?';
	
	if(projectionListXML != null){
		requestString += 'prlXML=' + projectionListXML + '&';
	}
	else if(projectionListDIR != null){
		requestString += 'prlDIR=' + projectionListDIR + '&';
	}
	else{
		alert('ProjectionList not specified!');
		return;
	}
	
	if(particleListXML != null){
		requestString += 'plXML=' + particleListXML + '&'
	}
	else{
		alert('No ParticleList specified!');
		return;
	}
	
	requestString += 'x=' + xCube + '&';
	requestString += 'y=' + yCube + '&';
	requestString += 'z=' + zCube + '&';
	
	requestString += 'sb=' + Before + '&';
	requestString += 'sa=' + After  + '&';
	
	requestString += 'xc=' + xCutPos + '&';
	requestString += 'yc=' + yCutPos + '&';
	requestString += 'zc=' + zCutPos + '&';
	
	requestString += 'ts=' + tomoBin;
	
	//display some additional boxes like - set filename of job, ...
	jobFile = prompt("Choose a Jobfile name. Provide full path! End with .sh!");

    if(jobFile.length > 0){
		requestString += '&jobFile=' + jobFile;
	} 
	
	//send request to server
	request(requestString,true);
	
	//alert(returnValue);
	
}


function createTomogramRecJobFromParameters(document){
	
	updateReconstructionCookiesFromPage();
	
	var projectionListDIR		= getCookie('pytom_reconstruction_projectionListDIR');
	var projectionListXML		= getCookie('pytom_reconstruction_projectionListXML');
	var tomogramName			= getCookie('pytom_reconstruction_tomogramName');
	var xCube 					= getCookie('pytom_reconstruction_xCube');
	var yCube 					= getCookie('pytom_reconstruction_yCube');
	var zCube 					= getCookie('pytom_reconstruction_zCube');
	var binBefore				= getCookie('pytom_reconstruction_binBefore');
	var binAfter				= getCookie('pytom_reconstruction_binAfter');
	
	requestString = 'createReconstructionJob.py?';
	
	if(projectionListXML != null){
		requestString += 'prlXML=' + projectionListXML + '&';
	}
	else if(projectionListDIR != null){
		requestString += 'prlDIR=' + projectionListDIR + '&';
	}
	else{
		alert('ProjectionList not specified!');
		return;
	}
	
	if(tomogramName != null){
		requestString += 'tomo=' + tomogramName + '&'
	}
	else{
		alert('No tomogram file specified!');
		return;
	}
	
	requestString += 'x=' + xCube + '&';
	requestString += 'y=' + yCube + '&';
	requestString += 'z=' + zCube + '&';
	
	requestString += 'sb=' + binBefore + '&';
	requestString += 'sa=' + binAfter  ;
	
	
	//display some additional boxes like - set filename of job, ...
	jobFile = prompt("Choose a Jobfile name. Provide full path! End with .sh");

    if(jobFile.length > 0){
		requestString += '&jobFile=' + jobFile;
	} 
	
	//send request to server
	request(requestString,true);
	
	//alert(returnValue);
	
}
