

raw_read_template <- function(tmpl_loc)
{
	raw_string <- readChar(tmpl_loc, file.info(tmpl_loc)$size)
	return(raw_string)
}

raw_write_page <- function(targetfile, str_content)
{
	cat(str_content, file=targetfile, sep='', fill=FALSE, append=FALSE)
	return(NULL)
}

add_new_content_page <- function(sourcefile, targetfolder, targetname)
{
	targetfile <- paste(targetfolder, targetname, sep='/')
	file.copy(sourcefile, targetfile, overwrite=TRUE)
	return(targetfile)
}

replace_css_menu_top_entry <- function(entry_name, color, css_styles)
{
	menu_top_entry <- 'nav ul li#<replaceName>:target a, nav ul li#<replaceName>:target > ul li a { background-color: #<replaceColor>;}'
	new_menu_entry <- gsub('<replaceName>', entry_name, x=menu_top_entry, fixed=TRUE)
	new_menu_entry <- gsub('<replaceColor>', color, x=new_menu_entry, fixed=TRUE)
	css_replacement <- paste(new_menu_entry, '/*replaceMenu*/', sep='\n\n')
	css_styles <- gsub('/*replaceMenu*/', css_replacement, x=css_styles, fixed=TRUE)
	return(css_styles)
}

replace_toplevel_menu_nosub <- function(entry_id, entry_html_ref, entry_text, nav_tmpl)
{
	replacement <- paste('<li id="', entry_id, '"><a href="', entry_html_ref,'">', entry_text, '</a></li>\n', sep='')
	replacement <- paste(replacement, '<!--replaceTopLevelMenu-->', sep='')
	new_nav <- gsub('<!--replaceTopLevelMenu-->', replacement, x=nav_tmpl, fixed=TRUE)
	return(new_nav)
}

replace_toplevel_menu_sub <- function(entry_id, entry_html_ref, entry_text, nav_tmpl)
{
	sub_id <- paste('<!--replaceSub', entry_id, '-->', sep='')
	replacement <- paste('<li id="', entry_id, '"><a href="#', entry_html_ref,'" target="navigation">', entry_text, '</a>\n<ul>\n', sub_id,'\n</ul>\n</li>\n', sep='')
	replacement <- paste(replacement, '<!--replaceTopLevelMenu-->', sep='')
	new_nav <- gsub('<!--replaceTopLevelMenu-->', replacement, x=nav_tmpl, fixed=TRUE)
	return(list(navstr=new_nav, subid=sub_id))
}

replace_submenu <- function(entry_id, entry_html_ref, entry_text, replace_id, extend, nav_tmpl)
{
	replacement <- paste('<li id="', entry_id, '"><a href="', entry_html_ref,'">', entry_text, '</a></li>\n', sep='')
	if (extend)
	{
		replacement <- paste(replacement, replace_id, sep='\n')
	}
	new_nav <- gsub(replace_id, replacement, x=nav_tmpl, fixed=TRUE)
	return(new_nav)
}

create_content_paragraph <- function(para_text)
{
	paragraph <- paste('<p>', para_text, '</p>', sep='')
	return(paragraph)
}

add_content_hrline <- function(content_page)
{
	replacement <- paste('<hr>', '<!--replaceContentElement-->', sep='\n')
	new_content_page <- gsub('<!--replaceContentElement-->', replacement, x=content_page, fixed=TRUE)
	return(new_content_page)
}

add_formatted_content <- function(form_cont, content_page)
{
	replacement <- paste(form_cont, '<!--replaceContentElement-->', sep='\n\n')
	new_content_page <- gsub('<!--replaceContentElement-->', replacement, x=content_page, fixed=TRUE)
	return(new_content_page)
}

create_table <- function(table_caption, table_header, table_data)
{
	stopifnot(length(table_header) == ncol(table_data))
	final_table <- paste('<table border="1" rules="groups" cellspacing="3" cellpadding="3">\n<caption>', table_caption, '</caption>\n<thead align="center" valign="middle">\n<tr>\n', sep='')
	n <- length(table_header)
	# the following is sometimes necessary since - apparently - labels are read as levels/categories in R
	# and converted to a digit when as.character is called... JFC!
	table_header <- as.vector(unlist(Map(toString, table_header)))
	header_part <- paste(c(rep('<th>', n)), table_header, c(rep('</th>', n)), sep='', collapse='\n')
	final_table <- paste(final_table, header_part, '\n', '</tr>\n</thead>\n<tbody align="center" valign="middle">', sep='')
	for (i in 1:nrow(table_data))
	{
		this_row <- as.vector(unlist(Map(toString, table_data[i,])))
		one_row <- paste(c(rep('<td>', n)), this_row, c(rep('</td>', n)), sep='', collapse='\n')
		final_table <- paste(final_table, '<tr>', one_row, '</tr>', sep='\n')
	}
	final_table <- paste(final_table, '</tbody>', '</table>', sep='\n')
	return(final_table)
}

create_list <- function(list_entries)
{
	final_list <- '<ul>'
	for (i in 1:length(list_entries))
	{
		one_entry <- paste('<li>', list_entries[i], '</li>')
		final_list <- paste(final_list, one_entry, sep='\n')
	}
	final_list <- paste(final_list, '</ul>', '\n', sep='\n')
	return(final_list)
}

wrap_section <- function(sectionbody, section_header, header_type)
{
	header_open <- list('<h1>', '<h2>', '<h3>', '<h4>')
	header_close <- list('</h1>', '</h2>', '</h3>', '</h4>')
	header_types <- c('h1', 'h2', 'h3', 'h4')
	stopifnot(header_type %in% header_types)
	names(header_open) <- header_types
	names(header_close) <- header_types
	sec_header <- paste(header_open[[header_type]], section_header, header_close[[header_type]], sep='')
	wrapped_section <- paste('<section>', sec_header, sectionbody, '</section>', sep='\n')
	return(wrapped_section)
}

wrap_visibility_toggle <- function(divbody, divid, linktext)
{
	anchor <- '<a href="#" onclick="toggle_visibility(\'<replaceID>\');">Click to show/hide: <replaceLinkText></a><br>'
	div_open <- '<div id="<replaceID>" style="display:block">'
	anchor <- gsub('<replaceID>', divid, x=anchor, fixed=TRUE)
	anchor <- gsub('<replaceLinkText>', linktext, x=anchor, fixed=TRUE)
	div_open <- gsub('<replaceID>', divid, x=div_open, fixed=TRUE)
	wrapped_vis <- paste(anchor, div_open, divbody, '</div>', sep='\n')
	return(wrapped_vis)
}
