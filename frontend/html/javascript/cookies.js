/**
 * 
 */



function setCookie(cookieName,value,erase){

	if (!erase) {
		//create cookie for a day
		var date = new Date();
		date.setTime(date.getTime()+(1*24*60*60*1000));
		var expires = "; expires="+date.toGMTString();
	}
	else{
		//erase cookie
		var date = new Date();
		date.setTime(date.getTime()-(1*24*60*60*1000));
		var expires = "; expires="+date.toGMTString();
	}

	document.cookie = cookieName+"="+value+expires+"; path=/";
	
}


function getCookie(cookieName) {
	var cookieNameEQ = cookieName + "=";
	var values = document.cookie.split(';');
	
	for(var i=0;i < values.length;i++) {
		var value = values[i];
		
		while (value.charAt(0)==' ') 
			value = value.substring(1,value.length);
		
		if (value.indexOf(cookieNameEQ) == 0) 
			return value.substring(cookieNameEQ.length,value.length);
		
	}
	return null;
}

function eraseCookie(cookieName) {
	createCookie(cookieName,"",True);
}