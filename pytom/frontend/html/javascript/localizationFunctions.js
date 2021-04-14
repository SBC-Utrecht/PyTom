function updateLocalizationCookiesFromPage(){
	var angleList = document.getElementById('AngleList');
	var globalSamplingFile = angleList[angleList.selectedIndex].value;
	setCookie('pytom_localization_angleFile',globalSamplingFile,false);
	
	var destination = document.getElementById('Destination').value;
	setCookie('pytom_localization_destination',destination,false);
	
	var lowestFrequency 	= document.getElementById('LowestFrequency').value;
	var highestFrequency 	= document.getElementById('HighestFrequency').value;
	var smooth				= document.getElementById('FilterSmooth').value;
	
	setCookie('pytom_localization_lowestFrequency',lowestFrequency);
	setCookie('pytom_localization_highestFrequency',highestFrequency);
	setCookie('pytom_localization_filterSmooth',smooth);
	
	var wedgeAngle1	= document.getElementById('Wedge1').value;
	var wedgeAngle2	= document.getElementById('Wedge2').value;
	
	setCookie('pytom_localization_wedgeAngle1',wedgeAngle1);
	setCookie('pytom_localization_wedgeAngle2',wedgeAngle2);
	
	var xCubes = document.getElementById('xCubes').value;
	var yCubes = document.getElementById('yCubes').value;
	var zCubes = document.getElementById('zCubes').value;
	
	setCookie('pytom_localization_xCubes',xCubes);
	setCookie('pytom_localization_yCubes',yCubes);
	setCookie('pytom_localization_zCubes',zCubes);
	
}

function createLocalizationJobFromParameters(){
	alert('Implement createLocalizationJobFromParameters');
	
}

function eraseLocalizationCookies(){
	setCookie('pytom_localization_tomogramPath','',true);
	setCookie('pytom_localization_reference','',true);
	setCookie('pytom_localization_mask','',true);
	setCookie('pytom_localization_angleFile','',true);
	setCookie('pytom_localization_destination','',true);
	setCookie('pytom_localization_lowestFrequency','',true);
	setCookie('pytom_localization_highestFrequency','',true);
	setCookie('pytom_localization_smooth','',true);
	setCookie('pytom_localization_wedgeAngle1','',true);
	setCookie('pytom_localization_wedgeAngle2','',true);
	setCookie('pytom_localization_xCubes','',true);
	setCookie('pytom_localization_yCubes','',true);
	setCookie('pytom_localization_zCubes','',true);
}

function resetLocalizationValues(){
	alert('Implement resetLocalizationValues');
	
}

function checkLocalizationParameters(theOverlay){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
	
		var tomogramPath = getCookie('pytom_localization_tomogramPath');
		var referencePath = getCookie('pytom_localization_reference');
		var maskPath = getCookie('pytom_localization_mask');
		var angleFile = getCookie('pytom_localization_angleFile');
		var destination = getCookie('pytom_localization_destination');
		var lowestFrequency = getCookie('pytom_localization_lowestFrequency');
		var highestFrequency = getCookie('pytom_localization_highestFrequency');
		var filterSmooth = getCookie('pytom_localization_filterSmooth');
		var wedgeAngle1 = getCookie('pytom_localization_wedgeAngle1');
		var wedgeAngle2 = getCookie('pytom_localization_wedgeAngle2');
		var xCubes = getCookie('pytom_localization_xCubes');
		var yCubes = getCookie('pytom_localization_yCubes');
		var zCubes = getCookie('pytom_localization_zCubes');
		
		theOverlay.style.display = "";
	
		var emptyLine 			= "<tr><td align='center' align='middle'><font class='checkTextName'></font></td>"+
		  "<td><font class='checkTextValue'>&nbsp;</font></td>" +
		  "</tr>";
	
		content = 	"<center>"+
					"<table>" +
					"<tr><td align='center' align='middle'><font class='checkTextName'>Tomogram</font></td>"+
					"<td><font class='checkTextValue'>&nbsp;" + tomogramPath + "</font></td>" +
					"</tr>";
					
		content +=  emptyLine;
		
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Reference</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + referencePath + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
	
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Mask</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + maskPath + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Angle file</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + angleFile + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Destination</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + destination + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Lowest Frequency</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + lowestFrequency + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Highest Frequency</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + highestFrequency + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Filter Smoothing</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + filterSmooth + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Wedge Angle 1</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + wedgeAngle1 + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Wedge Angle 2</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + wedgeAngle2 + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>X Cubes</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + xCubes + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Y Cubes</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + yCubes + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		content += 	"<center>"+
		"<table>" +
		"<tr><td align='center' align='middle'><font class='checkTextName'>Z Cubes</font></td>"+
		"<td><font class='checkTextValue'>&nbsp;" + zCubes + "</font></td>" +
		"</tr>";
		
		content +=  emptyLine;
		
		content += 	"<tr><td></td><td><input type='button' value='Cancel' onClick='return hideOverlay()'/></td></tr></table></center>";
		
		theOverlay.innerHTML = content; 
	}else{
		hideOverlay();
	}
return false;
	
}

function createLocalizationJobFromParameters(document){

    var tomogramPath = getCookie('pytom_localization_tomogramPath');
    var referencePath = getCookie('pytom_localization_reference');
    var maskPath = getCookie('pytom_localization_mask');
    var angleFile = getCookie('pytom_localization_angleFile');
    var destination = getCookie('pytom_localization_destination');
    var lowestFrequency = getCookie('pytom_localization_lowestFrequency');
    var highestFrequency = getCookie('pytom_localization_highestFrequency');
    var filterSmooth = getCookie('pytom_localization_filterSmooth');
    var wedgeAngle1 = getCookie('pytom_localization_wedgeAngle1');
    var wedgeAngle2 = getCookie('pytom_localization_wedgeAngle2');
    var xCubes = getCookie('pytom_localization_xCubes');
    var yCubes = getCookie('pytom_localization_yCubes');
    var zCubes = getCookie('pytom_localization_zCubes');

    var requestString = 'createLocalizationJob.py?';
    requestString += 'tomo=' + tomogramPath + '&';
    requestString += 'ref=' + referencePath + '&';
	requestString += 'mask=' + maskPath + '&';
    requestString += 'angle=' + angleFile + '&';
    requestString += 'dest=' + destination + '&';
    requestString += 'low=' + lowestFrequency + '&';
    requestString += 'high=' + highestFrequency + '&';
    requestString += 'smooth=' + filterSmooth + '&';
    requestString += 'w1=' + wedgeAngle1 + '&';
    requestString += 'w2=' + wedgeAngle2 + '&';
    requestString += 'x=' + xCubes + '&';
    requestString += 'y=' + yCubes + '&';
    requestString += 'z=' + zCubes;

    //display some additional boxes like - set filename of job, ...
	jobFile = prompt("Choose a Jobfile name. Provide full path! End with .xml");

    if(jobFile.length > 0){
		requestString += '&jobFile=' + jobFile;
	} 
	
	//send request to server
	request(requestString,true);
	
	//alert(returnValue);

}
