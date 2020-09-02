SELECT  count(distinct p.id) FROM sndg5.bioresources_person p , sndg5.bioresources_affiliation a ,
	sndg5.bioresources_affiliation_organizations ao , sndg5.bioresources_organization o,
    sndg5.bioresources_resourcerelation r
WHERE p.id = a.author_id AND ao.affiliation_id = a.id AND o.id = ao.organization_id
	AND o.country = "Argentina" AND a.publication_id = r.source_id AND LOWER(r.role) LIKE  'publication_%';

SELECT role,count(distinct target_id) FROM sndg5.bioresources_resourcerelation
group by role;

SELECT v.value,count(*) as cant FROM sndg5.bioresources_resourcepropertyvalue v, sndg5.bioresources_resourceproperty p
WHERE (p.term_id = 30 or p.term_id = 96) AND p.id = v.property_id
GROUP BY v.value
ORDER BY cant DESC;

SELECT count(*) as cant FROM sndg5.bioresources_resourcepropertyvalue v, sndg5.bioresources_resourceproperty p
WHERE p.term_id in (18,19,20) AND p.id = v.property_id;

SELECT t.identifier,t.term_id, count(*) as cant
FROM sndg5.bioresources_resourcepropertyvalue v,
 sndg5.bioresources_resourceproperty p,
 sndg5.term t
WHERE v.property_id = p.id AND t.term_id = p.term_id
GROUP BY p.term_id
ORDER BY cant DESC;