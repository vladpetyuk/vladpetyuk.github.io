---
layout: page
title: posts
subtitle: a place to vent
---

<!--## Posty Postyarychiii-->

<!--
{% for post in site.posts %}
{% if post.categories contains 'category1' %}
  <div class="post">
    <h3 class="title"><a href="{{ post.url }}">{{ post.title }}</a></h3>
    <p class="meta">{{ post.date | date_to_string}}</p>
    <div class="entry">
     {{ post.content | strip_html | truncatewords: 100 }}
    </div>
 </div>
{% endif %}
{% endfor %}
-->


{% for post in site.posts %}
{% if post.categories contains 'category1' %}
<div class="posts-list">

  <article class="post-preview">
    <a href="{{ post.url | prepend: site.baseurl }}">
	  <h2 class="post-title">{{ post.title }}</h2>
	
	  {% if post.subtitle %}
	  <h3 class="post-subtitle">
	    {{ post.subtitle }}
	  </h3>
	  {% endif %}  
    </a>

    <p class="post-meta">
      Posted on {{ post.date | date: "%B %-d, %Y" }}
    </p>
  
    <div class="post-entry">
      {{ post.content | truncatewords: 50 | strip_html}}
	  <a href="{{ post.url | prepend: site.baseurl }}" class="post-read-more">[Read&nbsp;More]</a>
    </div>
  
   </article>
</div>

{% if paginator.total_pages > 1 %}
<ul class="pager main-pager">
  {% if paginator.previous_page %}
  <li class="previous">
    <a href="{{ paginator.previous_page_path | prepend: site.baseurl | replace: '//', '/' }}">&larr; Newer Posts</a>
  </li>
  {% endif %}
  {% if paginator.next_page %}
  <li class="next">
    <a href="{{ paginator.next_page_path | prepend: site.baseurl | replace: '//', '/' }}">Older Posts &rarr;</a>
  </li>
  {% endif %}
</ul>
{% endif %}

{% endif %}
{% endfor %}





<!--
{% for post in site.posts %}
{% if post.categories contains 'category1' %}
   * {{ post.date | date_to_string }} &raquo; [ {{ post.title }} ]({{ post.url }})
{% endif %}
{% endfor %}
-->

<!--
layout: default
title: Blog1
category: blog1
-->
