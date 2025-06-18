mail='18706726158@163.com';%自己邮箱地址    需要自己设置
password='JHURFRQFMKNJPNVK';%%是163邮箱开通SMTP时出现的密码    需要自己设置

%%服务器设置
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.163.com');%163的smtp应该是固定的
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props=java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','456');

%%发邮件
accuracy=0.888;
receiver='550399752@qq.com';%收件人邮箱
mailtitle='pattern recognition';%邮件主题
mailcontent=['mission compeleted!',...    %邮件内容
    'accuracy=',num2str(accuracy)];       %邮件内容
sendmail(receiver,mailtitle,mailcontent);   %发送