

    function loadFile() {
        var input, file, fr, result;
        var x_length, y_length, z_length;



        if (typeof window.FileReader !== 'function') {
            bodyAppend("p", "The file API isn't supported on this browser yet.");
            return;
        }

        input = document.getElementById('fileinput');
        if (!input) {
            bodyAppend("p", "Um, couldn't find the fileinput element.");
        }
        else if (!input.files) {
            bodyAppend("p", "This browser doesn't seem to support the `files` property of file inputs.");
        }
        else if (!input.files[0]) {
            bodyAppend("p", "Please select a file before clicking 'Load'");
        }
        else {
            file = input.files[0];
            fr = new FileReader();
            //fr.onload = receivedText;
            //fr.readAsText(file);

            fr.onload = receivedBinary;
            fr.readAsBinaryString(file);
        }
        
/*        function receivedText() {
            showResult(fr, "Text");

            fr = new FileReader();
            fr.onload = receivedBinary;
            fr.readAsBinaryString(file);
        }
*/      
 
        function receivedBinary() {
            showResult(fr, "Binary");
            sizeImage(fr, "Binary", x_length, y_length, z_length);
            //showSlid(fr, "Binary",x_length, y_length, z_length);
            floatArray = convert2Float(fr,'Binary');
            floatArray = normalize(floatArray);
        }
        
        return floatArray;
    }
    
    function normalize(floatArray) {
        var length = floatArray.length;
        var max, min;
        // find the max and min value
        for (var i = 0; i < length; ++i) {
            if (i == 0) {
                max = floatArray[0];
                min = floatArray[0];
            }
            if (max > floatArray[i])
                max = floatArray[i];
            if (min < floatArray[i])
                min = floatArray[i];
        }
        
        // normalization
        for (var i = 0; i < length; ++i) {
            floatArray[i] = (floatArray[i] - min)/(max-min);
        }
        return floatArray;
    }

    function showResult(fr, label) {
        var markup, result, n, aByte, byteStr;

        markup = [];
        result = fr.result;
        // for (n = 0; n < result.length; ++n) {
        for (n = 0; n < 100; ++n) {
            aByte = result.charCodeAt(n);
            byteStr = aByte.toString(10);
            if (byteStr.length < 2) {
                byteStr = "0" + byteStr;
            }
            markup.push(byteStr);
        }
        bodyAppend("p", label + " (" + result.length + "):");
        bodyAppend("pre", markup.join(" "));
    }
    
    function convert2Float(fr,label) {
        var markup, result, n, aByte, byteStr;

        markup = [];
        result = fr.result;
        // for (n = 0; n < result.length; ++n) {
        var floatArray = new Array();
        var head = 512;
        var typeLength = 4;
        var result = fr.result;
        var a, b, c, d, value;
        for (n = 0; n < result.length; n= n + typeLength) {
            a = result.charCodeAt(n+head);
            b = result.charCodeAt(n+1+head);
            c = result.charCodeAt(n+2+head);
            d = result.charCodeAt(n+3+head);
            value = a | b<<8 | c << 16| d << 24;
            value = hex2float(value);
            floatArray.push(value);
            
          //  valueString = value.toString();
          //  markup.push(valueString);           
        }
        //bodyAppend("p", label + " (" + result.length + "):");
        //bodyAppend("pre", markup.join(" "));
        return floatArray;
    }
    
    function hex2float(num) {
        var sign = (num & 0x80000000) ? -1 : 1;
        var exponent = ((num >> 23) & 0xff) - 127;
        var mantissa = 1 + ((num & 0x7fffff) / 0x7fffff);
        return sign * mantissa * Math.pow(2, exponent);
    }

    
    function middleSlice(floatArray, x_length, y_length, z_length) {
        var slice = new Array();
        var offset, n, i = 0;
        var zAix = 1/2 * z_length;
        offset = x_length * y_length * zAix;

        for (n = 0; n < x_length * y_length; n++ ) {
            //slice[i] = dataArray[n + offset];
            slice[i] = Math.random();
            i++;
        }
        return slice;

        //return dataArray.slice(offset,x_length * y_length);
    }
   

    function randomSlice(x_length, y_length, z_length) {
        var slice = new Array();
        var offset, n, i = 0;
        var zAix = 1/2 * z_length;
        offset = x_length * y_length * zAix;

        for (n = 0; n < x_length * y_length; n++ ) {
            slice[i] = Math.random();
            i++;
        }
        return slice;
    }



    function sizeImage(fr, label, x_length, y_length, z_length) {
        var markup,result, n, byteStr;
        var a, b, c ,d;
        markup = [];
        result = fr.result;
        for (n = 4; n < 16; n = n + 4) {
            a = result.charCodeAt(n);
            b = result.charCodeAt(n+1);
            c = result.charCodeAt(n+2);
            d = result.charCodeAt(n+3);
            if (n == 4) {
                 x_length = a | b<<8 | c << 16| d << 24;
                 byteStr = x_length.toString(10);
                 markup.push(byteStr);
            }
            if (n == 8) {
                y_length = a | b<<8 | c << 16| d << 24;
                byteStr = y_length.toString(10);
                markup.push(byteStr);
            }
            if ( n == 12) {
                z_length = a | b<<8 | c << 16| d << 24;
                byteStr = z_length.toString(10);
                markup.push(byteStr);
            }
    
        }
        bodyAppend("p", label + " (" + result.length + "):");
        bodyAppend("pre", markup.join(" "));
    }
/*
    function showSlid(fr, label, x_length, y_length, z_length) {
        var headSize = 4+12+80+40;
        var slidSize = x_length * y_length;
        var markup, result, n, aByte, byteStr;
        markup = [];
        result = fr.result;
        for (n = headSize; n < headSize+1000; ++n){
            aByte = result.charCodeAt(n);
            byteStr = aByte.toString(10);
            if (byteStr.length < 2) {
                byteStr = "0" + byteStr;
            }
            markup.push(byteStr);
        }
        bodyAppend("p", label + " (" + result.length + "):");
        bodyAppend("pre", markup.join(" "));
    }
*/



    function bodyAppend(tagName, innerHTML) {
        var elm;

        elm = document.createElement(tagName);
        elm.innerHTML = innerHTML;
        document.body.appendChild(elm);
    }
    
    