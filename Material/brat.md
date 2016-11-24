# 安装
## 在ubuntu安装apche2及支持cgi
1. 安装apche2
```
sudo apt-get install apache2
```
2. 安装CGI模块
```
sudo vim /etc/apache2/mods-enabled/mime.load
```
添加
```
LoadModule cgi_module /usr/lib/apache2/modules/mod_cgi.so
```
3. 修改虚拟机
```
sudo vim /etc/apache2/sites-available/000-default.conf
```
添加
```
        <Directory "/var/www/html">
                AllowOverride Options Indexes FileInfo Limit
                Options +ExecCGI -MultiViews +SymLinksIfOwnerMatch
                AddType application/xhtml+xml .xhtml
                AddType font/ttf .ttf
                        # For CGI support
                AddHandler cgi-script .cgi
                        # Comment out the line above and uncomment the line below for FastCGI
                        #AddHandler fastcgi-script fcgi
        </Directory>
```
取消注释
```
#Include conf-available/serve-cgi-bin.conf
```
