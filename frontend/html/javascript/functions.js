var xmlhttp = null;

function request(requestString,isAsynchronous){
	
	if(xmlhttp == null){
		if (window.XMLHttpRequest)
		{// code for IE7+, Firefox, Chrome, Opera, Safari
			xmlhttp=new XMLHttpRequest();
		}
		else
		{// code for IE6, IE5
			xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
		}
	}

	try
	{
		if(isAsynchronous && xmlhttp.onreadystatechange == null){
			xmlhttp.overrideMimeType('text/xml');
			xmlhttp.onreadystatechange = handleAsynchronousRequests;
		}
		xmlhttp.open("GET",requestString,isAsynchronous);
		xmlhttp.send();
		
	}
	catch(e){
		alert(e);
	}
	
	if(isAsynchronous)
		return dispatchHTTPReply(xmlhttp.responseXML);
	else
		return xmlhttp.responseText;
	//return xmlhttp.responseXML;
}

function requestSynchronous(requestString){
	return request(requestString,false);
}


function loadPage(pageName){
	return request(pageName);
}

function handleAsynchronousRequests(){
	if (xmlhttp.readyState == 4){
		
		var responseXML = null;
		if(xmlhttp.responseXML != null){
			responseXML = xmlhttp.responseXML;
		}else{
			responseXML = stringToXML(xmlhttp.responseText);
		}
		dispatchHTTPReply(responseXML);
		//displayJobFileCreatedResponse(xmlhttp.responseText);
		
	}
}

function dispatchHTTPReply(responseXML){
	
	var file 	= responseXML.getElementsByTagName('File');
	var error 	= responseXML.getElementsByTagName('Error');
	
	if(file.length > 0){
		//parse file response
		var type 	=file[0].getAttribute('type');
		var path 	=file[0].getAttribute('path');
		var status	=file[0].getAttribute('status');
		
		if(type in {'FileExistance':'', 'DirExistance':''}){
			return;
		}
		else if(type in {'ReconstructionJob':'', 'MCOEXMXJob':'','MCOACJob':'', 'LocalizationJob':'','AlignmentJob':'',}){
			displayMessage("File " + type + " was " + status + "! Path: " + path);
		}
	}
	
	if(error.length > 0){
		//parse error response
		var message = error[0].getAttribute('message');
		displayMessage("Error: " + message);
	}
	
}

function displayMessage(message){
	theOverlay = document.getElementById('overlayBox');
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		theOverlay.innerHTML = "<table width='100%' height='100%'>" +
		   "<tr><td align='center' align='middle' width='100%' height='100%'>" +
		   "<div id='overlayText' class='overlayText'><font class='checkTextName' >"+ message +" </font></div></td></tr><tr><td align='center' align='middle' width='100%' height='100%'>"+ 
		   "<input type='button' value='Cancel' onClick='hideOverlay()'/>" +
		   "</td></tr></table>";
	}
}


function waitMessageRunFunction(message,functionPointer){
	
	functionPointer(request(message));
	
}

function setHTMLFromPage(message,target){
	
	//alert(message.substring(0,100));
	target.innerHTML = request(message,false);
	//alert(target.innerHTML.substring(0,100));
}

function setImageFromPage(message,target){
	//alert(target);
	target.src = request(message,false);
	//alert(target.innerHTML.substring(0,100));
}

function dialogOverlay(theOverlay,message){
	
	setCookie('temporaryDialogCookie','',false);
	
	setCookieFromOverlay(theOverlay,message,'temporaryDialogCookie');
	
	value = '';
	
	while(value == ''){
		value = getCookie('temporaryDialogCookie');
	}
	
	setCookie('temporaryDialogCookie','',true);
	
	return value;
	
}

function showOverlayFileUpload(theOverlay,serverPage,asParameter){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		theOverlay.innerHTML = "<table width='100%' height='100%'>" +
							   "<tr><td align='center' valign='middle' width='100%' height='100%'>" +
							   "<form><input type='file' multiple='multiple' id='fileUpload'/>" +
							   "<div class='smallButton' onclick='alert('asdf')'>Load</div>&nbsp;" +
							   "<input type='button' value='Cancel' onClick='return hideOverlay()'/></form>" +
							   "</td></tr></table>";
	}else{
		hideOverlay();
	}
	return false;
}



function showOverlayImage(theOverlay,imageName){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		theOverlay.innerHTML = "<table width='100%' height='100%'>" +
							   "<tr><td align='center' valign='middle' width='100%' height='100%'>" +
							   "<img src='"+imageName+"'/>" +
							   "<input type='button' value='Cancel' onClick='return hideOverlay()'/>" +
							   "</td></tr></table>";
	}else{
		hideOverlay();
	}
	return false;
}


function setCookieFromOverlay(theOverlay,message,cookieName){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		theOverlay.innerHTML = "<table width='100%' height='100%'>" +
							   "<tr><td align='center' align='middle' width='100%' height='100%'>" +
							   "<div id='overlayText' class='overlayText'>"+ message +"</div>"+ 
							   "<form><input type='text' id='inputText' size='100' maxlength='500'/>" +
					
							   "<br/><input type='button' value='Load' onClick='setCookie(\""+cookieName+"\",document.getElementById(\"inputText\").value,false);hideOverlay();'/>" +
							   
							   "<input type='button' value='Cancel' onClick='hideOverlay()'/></form>" +
							   "</td></tr></table>";
	}else{
		hideOverlay();
	}
	return false;
}


function imageFromFileInputOverlay(theOverlay,serverPage,asParameter,imageTargetName,message,cookieName){
	//show overlay
	//request content when clicked, display image 
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		
		theOverlay.innerHTML = "<table width='100%' height='100%' style='opacity: 1.0;'>" +
							   "<tr><td align='center' align='middle' width='100%' height='100%'>" +
							   "<div id='overlayText' class='overlayText' style='opacity: 1.0;'>"+ message +"</div>"+ 
							   "<form><input type='text' id='filePath' size='100' maxlength='500'/>" +

							   "<br/><input type='button' value='Load' style='opacity: 1.0;' onClick='setImageFromPage(\""+ serverPage+"\" + \"?\" + \""+asParameter+"\" + \"=\"  + document.getElementById(\"filePath\").value , document.getElementById(\""+imageTargetName+"\") );setCookie(\""+cookieName+"\",document.getElementById(\"filePath\").value,false);hideOverlay();'/>" +
							   
							   "<input type='button' value='Cancel' style='opacity: 1.0;' onClick='hideOverlay()'/></form>" +
							   "</td></tr></table>";
		
	}else{
		hideOverlay();
	}
	return false;
}


function innerHTMLFromFileInputOverlay(theOverlay,serverPage,asParameter,htmlTargetName,message,cookieName){
	//request content when clicked, fill destination with content
	/**
	* Displays overlay, sends a request to a server page with some parameters, fills an html object with content, sets a cookie to some value 
	* @param {DomObject} theOverlay 
	* @return {String}   Returns a string value containing name and greeting
	*/
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		
		theOverlay.innerHTML = "<table width='100%' height='100%' style='opacity: 1.0;'>" +
							   "<tr><td align='center' align='middle' width='100%' height='100%'>" +
							   "<div id='overlayText' class='overlayText' style='opacity: 1.0;'>"+ message +"</div>"+ 
							   "<form><input type='text' id='filePath' size='100' maxlength='500'/ style='opacity: 1.0;'>" +

							   "<br/><input type='button' value='Load' style='opacity: 1.0;' onClick='setHTMLFromPage(\""+ serverPage+"\" + \"?\" + \""+asParameter+"\" + \"=\"  + document.getElementById(\"filePath\").value , document.getElementById(\""+htmlTargetName+"\") );setCookie(\""+cookieName+"\",document.getElementById(\"filePath\").value,false);hideOverlay();'/>" +
							   
							   "<input type='button' value='Cancel' style='opacity: 1.0;' onClick='hideOverlay()'/></form>" +
							   "</td></tr></table>";
		
	}else{
		hideOverlay();
	}
	return false;
}

function overlayWaitMessageRunFunction(theOverlay,message,serverPage,asParameter,functionPointer){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		
		theOverlay.style.display = "";
		
		theOverlay.innerHTML = "<table width='100%' height='100%' style='opacity: 1.0;'>" +
							   "<tr><td align='center' align='middle' width='100%' height='100%'>" +
							   "<div id='overlayText' class='overlayText' style='opacity: 1.0;'>"+ message +"</div>"+ 
							   "<form><input type='text' id='filePath' size='100' maxlength='500'/ style='opacity: 1.0;'>" +

							   "<br/><input type='button' value='Load' style='opacity: 1.0;' onClick='waitMessageRunFunction(\""+ serverPage+"\" + \"?\" + \""+asParameter+"\" + \"=\"  + document.getElementById(\"filePath\").value , \""+ functionPointer+ "\");hideOverlay();'/>" +
							   
							   "<input type='button' value='Cancel' style='opacity: 1.0;' onClick='hideOverlay()'/></form>" +
							   "</td></tr></table>";
		
	}else{
		hideOverlay();
	}
	return false;
}
function showOverlayMessage(theOverlay,message){
	
	if(!theOverlay){
		hideOverlay();
	}
	else if(theOverlay.style.display == "none"){
		theOverlay.style.display = "";
		theOverlay.innerHTML = "<table width='100%' height='100%' style='opacity: 1.0;'><tr>" +
							   "<td align='center' valign='middle' width='100%' height='100%' style='opacity: 1.0;'>"+message+"<br><br>" +
							   	"<input type='button' style='opacity: 1.0;' value='Cancel' onClick='return hideOverlay()'/>" +
							   	"</td></tr></table>";
	}else{
		hideOverlay();
	}
	return false;
}

function hideOverlay(){
	var theOverlay = document.getElementById('overlayBox');
	theOverlay.style.display = "none";
	theOverlay.innerHTML = '';
	return false;
}

function getAPIDocumentationURL(){
	
	apiPath = loadPage('apiDocumentationURL.py');
	
	return apiPath;
	
}

function stringToXML(text){
    if (window.ActiveXObject){
      var newDoc=new ActiveXObject('Microsoft.XMLDOM');
      newDoc.async='false';
      newDoc.loadXML(text);
    } else {
      var parser=new DOMParser();
      var newDoc=parser.parseFromString(text,'text/xml');
    }
    return newDoc;
}

/*
function setCookie(cookieName,value){
	var c_value=escape(value) +  "expires="+);
	document.cookie= cookieName + "=" + value;
	
	
}
*/